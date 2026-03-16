// OLSRv2 集成功能测试：nhdp、rfc5444
// 路径：mysrc/olsr/test/test_function
// 编译：g++ -std=c++17 -I../include -o test_function test_function.cc ../src/olsrv2.cc ../src/olsrv2_state.cc ../src/olsrv2_tc.cc ../src/nhdp/nhdp_db.cc ../src/rfc5444/rfc5444_reader.cc ../src/rfc5444/rfc5444_writer.cc
// 运行：./test_function

#include <iostream>
#include <cassert>
#include "../include/olsrv2.h"
#include "../include/olsrv2_state.h"
#include "../include/nhdp/nhdp_db.h"
#include "../include/rfc5444/rfc5444_reader.h"
#include "../include/rfc5444/rfc5444_writer.h"

using namespace std;

void test_olsrv2_nhdp_integration() {
    // 初始化 OLSRv2 状态
    Olsrv2State state;
    // 添加邻居
    Ipv4Address neighbor("10.0.0.2");
    state.nhdpDb.addNeighbor(neighbor);
    assert(state.nhdpDb.isNeighbor(neighbor));
    // 触发邻居集变更，检查 OLSRv2 路由集是否响应
    state.updateRoutingSet();
    assert(state.routingSet.count(neighbor) > 0);
    cout << "[PASS] OLSRv2-NHDP 集成: 邻居同步与路由集更新" << endl;
}

void test_olsrv2_rfc5444_integration() {
    // 构造 OLSRv2 Hello 消息
    Olsrv2Hello hello;
    hello.originator = Ipv4Address("10.0.0.1");
    hello.neighbors.push_back(Ipv4Address("10.0.0.2"));
    // 编码为 RFC5444
    std::vector<uint8_t> buffer;
    assert(rfc5444::encodeHello(hello, buffer));
    // 解码回 OLSRv2 Hello
    Olsrv2Hello decoded;
    assert(rfc5444::decodeHello(buffer, decoded));
    assert(decoded.originator == hello.originator);
    assert(decoded.neighbors == hello.neighbors);
    cout << "[PASS] OLSRv2-RFC5444 集成: Hello 编解码" << endl;
}

int main() {
    test_olsrv2_nhdp_integration();
    test_olsrv2_rfc5444_integration();
    cout << "[ALL PASS] OLSRv2 功能集成测试完成" << endl;
    return 0;
}

