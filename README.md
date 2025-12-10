
## 更新日志（12.10）
- 修改原有DT相关计算逻辑，更加简洁易用，不使用缓存，直接每次重算（增量更新不知道为什么CGAL库容易崩）
- 控制包完善转发逻辑（之前的转发包处理放在类似processNeighborSetNotificationChunk命名函数中，实际上没有处理到），见MDTRouting::ensureRouteForDatagram中对控制包的特殊处理，控制包逻辑和协议对比改了一些错误实现（如应该发送自己的dt邻居集合而不是帮别人计算dt邻居）
- 注意标记no use的逻辑部分，应该是没用的
- 后续完善：丢包的处理，原论文中没有处理什么时候应该判断不可达丢包，导致实际可能成环路，就是发不过去


## 关于框架的解释
Omnet++特有几种文件，ned，ini，msg，分别是接入仿真的模块，整个仿真的参数设置，数据包定义（msg文件会被自动处理为.h和.cc文件），以及代码除了本身的路由逻辑，还有部分debugger的逻辑（把日志发出来然后自动处理），分布在RPC文件夹和routing文件夹的部分文件中。和路由本身相关的肯定还是都在routing文件夹中，其中大量使用了框架接口，所有框架代码都在inet-4.5.4/src/inet中。
## 代码核心架构
- 总控制：MDTRouting.h/cc
- 转发控制（选路与更新）：ForwardProtocol.h/cc
- 邻居管理：Neighbor.h/cc（很多三角剖分的计算放在NeighborDT.h/cc里面了（使用CGAL库），计算逻辑应该没什么大问题，有过初步测试，但可能还是有bug）
- 初始化和维护相关：JoinProtocol和MaintenanceProtocol（这里很多逻辑是为了兼容处理，节点离开网络和加入网络，但是我还没有设计当节点发现自己即将失去最后一个邻居的特殊处理，所以认为初始化完成后节点状态一直是在网络中）
## 对包发送的说明
包分为数据包和控制包
### 数据包
使用这几个钩子可以在不同阶段获取数据包
```c++
    //when a node is going to send a packet or receiving a packet(the packet just arrives network layer)
    Result ensureRouteForDatagram(Packet *datagram,int state);
    virtual Result datagramPreRoutingHook(Packet *datagram) override //刚刚收到
    virtual Result datagramForwardHook(Packet *datagram) override  //决定转发之后
    virtual Result datagramPostRoutingHook(Packet *datagram) override { return ACCEPT; }//即将离开本节点
    virtual Result datagramLocalInHook(Packet *datagram) override //即将本地接收（本节点是目的地）
    //local process the packet(other layer in charge of it)
    virtual Result datagramLocalOutHook(Packet *datagram) override  //从本地发送出去的包（本地出发）
```
ensureRouteForDatagram中对包进行处理，查找下一跳，修改路由表
### 控制包
参见MDTControlPackets.msg文件，其中定义了所有控制包，注意 @customize(true);标记表示该包的数据成员使用了msg没有直接支持的数据结构，这种数据包有专门的头文件，其中添加了额外的数据成员，在MDTRouting.h/cc的总控制中，使用void MDTRouting::sendPacketTo(Packet *packet, const L3Address& destAddr, int destPort,double delay)函数发送控制包（个别控制包不是，但都使用的是socket.sendTo接口），由于路由模块继承了public UdpSocket::ICallback,所以可以使用如下接口接收数据包：
```c++
    /* UDP callback interface */
    virtual void socketDataArrived(UdpSocket *socket, Packet *packet) override;
    virtual void socketErrorArrived(UdpSocket *socket, Indication *indication) override;
    virtual void socketClosed(UdpSocket *socket) override;
```
对控制包的分流处理在函数：void processPacket(Packet *pk);中，在分流前有一个接收功率判定逻辑，通信节点自己有一个最小接收功率，然后人工定义了一个更加窄的阈值，如果超过该阈值，认为是脆弱信号，在KA和KAN消息的处理中要根据是否是脆弱信号，进行分支处理。分流之后包会在不同组件处理（TODO：有一些控制包会转发，但是目前的逻辑是重新构造一个包再发出去，对overhead的评估必然导致偏高）具体控制包：
```java
class MDTControlPacket extends FieldsChunk
{
    MDTControlPacketType packetType = static_cast<MDTControlPacketType>(-1);       // packet type, stored in one byte
}

class KeepAlive extends MDTControlPacket {//心跳，携带自己可以听到哪些节点
	//L3Address sourceAddr;
	//simtime_t timestamp;
	@customize(true);
	double coords[];
	//int dim;
	bool attached;
}
class KeepAliveNotification extends MDTControlPacket {//心跳，只有自己坐标
	double coords[];
	//int dim;
	bool attached;
}
class JoinRequest extends MDTControlPacket{//加入申请
    double coords[];
    //Coordinate nodeCoord;
    int dim;
    L3Address destAddr;
    L3Address sourceAddr;
    L3Address predAddr;
}
class JoinReply extends MDTControlPacket{//加入回复，携带计算得到的请求方dt邻居
    @customize(true);
    L3Address sourceAddr;
    L3Address predAddr;
    L3Address destAddr;
}
class NeighborSetRequest extends MDTControlPacket{//邻居请求
    double coords[];
    //Coordinate nodeCoord;
    int dim;
    L3Address sourceAddr;
    L3Address destAddr;
    L3Address predAddr;
}
class NeighborSetReply extends MDTControlPacket{//邻居回复
    @customize(true);
    L3Address sourceAddr;
    L3Address predAddr;
    L3Address destAddr;
}
class NeighborSetNotification extends MDTControlPacket{//弃用，因为KeepAliveNotification完全替代了其功能
    double coords[];
    int dim;
    L3Address sourceAddr;
    L3Address predAddr;
    L3Address destAddr;
}
class CoordsDiscoverRequest extends MDTControlPacket{//加快获取坐标，如果本地没有目的地坐标，
                                                    // 向单跳邻居请求坐标（不扩散），实际上这个包可以有很多拓展（如果本地
                                                    // 信息过期也可以发，目前本地已知节点Known集合不会过期）
    //L3Address destAddr;
    //int num;
    L3Address sourceAddr;
    L3Address targetAddr;
    int hopcount;
}
class CoordsDiscoverReply extends MDTControlPacket{//获取坐标的回复
    @customize(true);
    L3Address destAddr;
    L3Address sourceAddr;
    L3Address targetAddr;
}
```
ps：还有一个控制overhead的方法，所有坐标都使用float，不使用double，听说在对overhead更加严格的场景还可以使用半精度浮点数
## 心跳信号和本地节点集合更新
本地有三个集合，物理邻居（单跳可达），dt邻居，和已知节点（known）集合，KeepAlive和KeepAliveNotification控制包频率较高，刷新物理邻居，每一次刷新，如果有物理邻居因为过期导致被删除，随即会触发路由更新，在MDTRouting.h/cc中的方法：void MDTRouting::deleteexpiredneighbors(std::vector<L3Address> nodes)，使用本地路由探查，尝试获取下一跳，如果贪心可以得到下一跳自然最好，但是如果不行，那么使用MDT的本地探查，很可能导致二环路（对于目的地d，a认为应该先去b，b认为应该先去a），这里应该增加的逻辑是如果建立的是mdt下一个目的地，最好使用控制包发一次，看能不能跑通，同理前述的ensureRouteForDatagram也应该对dt链路进行确认。dt链路确认的工作本来是由MaintenanceProtocol负责的，但是它是固定周期维护，如果可以处理为只要dt链路扰动就直接触发维护，应该效果更好。dt邻居计算单独周期维护（如果用增量更新本地邻居图应该可以实时维护）
## 对数据包附加额外字段
MdtRelayChunk，标记当前数据包处于贪心转发还是mdt链路发送（如果是后者要优先送到mdt邻居，才可以绕过local min），还有前一个节点predaddr，用来简单诊断二环路
## 路由条目属性
目前只使用了两个属性，目的地和下一跳节点，没有使用metric来衡量链路质量（应该加上），没有使用原论文中的四元组<src,pred,succ,dest>,虽然有MDTData.h/cc规范这个属性，但是缺乏维护，也没有使用，但是如果要更加好地维护链路，这个维护也有必要
## 初始化
JoinProtocol和MaintenanceProtocol是项目最开始写的代码，可能有点难读，它们是为了适应完全零启动，但是实际上我们可以让各个节点启动的时候就已经知道其他所有节点的位置，然后算好全局的dt图向其他节点各发一个包，确认可以收到和更新路由表，然后再启动所有节点开始移动和发包，才开始写的时候没有想这么多，只是从复现论文的角度全把他们写进去了，绕了一些弯路
## 目前实现的不足
- 为了尽快实现，我是直接对照框架已有的路由协议，使用它们使用的接口，而没有去研究具体和其他层如何对接。
- 目前实现对动态性的适应还是不足，维护协议的适配做的不够（而且KeepAlive拓展后持有自己听到的节点集合，相当于也该承担维护责任，我没有让它承担维护，通过KA发现的新dt邻居没有去维护）
- 可能还有bug，虽然比较严重的都修复了，但是bug总改不完


