#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../include/nhdp/nhdp.h"
#include "../include/olsrv2_routing.h"
#include "../include/rfc5444/rfc5444.h"
#include "../message/OLSRv2Packet_m.h"

#include "inet/common/packet/Packet.h"
#include "inet/networklayer/common/L3Address.h"

#include "test_runner.h"

using namespace mysrc::olsrv2;

namespace {
int g_failures = 0;

void expectTrue(bool cond, const std::string& msg)
{
    if (!cond) {
        ++g_failures;
        std::cerr << "[FAIL] " << msg << "\n";
    }
}

const RouteTuple *findRoute(const std::vector<RouteTuple>& routes, const MainAddress& dst)
{
    for (const auto& r : routes) {
        if (r.destination == dst)
            return &r;
    }
    return nullptr;
}

void testEndToEndHelloTcRoute()
{
    // 集成测试目标：
    // 1) 通过 NHDP 处理 HELLO，构建一跳对称邻居（N1<->N2）
    // 2) 通过 OLSRv2Core 从 NHDP 数据生成 TC（N1 广播其对称邻居 N2）
    // 3) 通过 OLSRv2Core 处理收到的 TC，在 S 节点生成到 N2 的路由（下一跳为 N1）
    // 4) 通过 RFC5444 的时间 TLV 编解码辅助，验证时间/有效期相关工具链可用（协议集成点之一）

    const MainAddress S("10.0.0.1");
    const MainAddress N1("10.0.0.2");
    const MainAddress N2("10.0.0.3");

    const double now = 100.0;

    // RFC5444 时间 TLV 编解码：模拟“把协议定时参数映射到报文编码/解析”的集成路径
    const uint64_t vtime_ms = 5000;
    const uint8_t encoded_vtime = mysrc::olsrv2::rfc5444::rfc5497_timetlv_encode(vtime_ms);
    const uint64_t decoded_vtime_ms = mysrc::olsrv2::rfc5444::rfc5497_timetlv_decode(encoded_vtime);
    expectTrue(decoded_vtime_ms > 0, "RFC5444 timetlv encode/decode must return positive value");

    // 节点 N1：运行 NHDP，接收来自 N2 的 HELLO，形成 N1<->N2 对称链路
    Nhdp nhdpN1;
    nhdpN1.setOriginator(inet::L3Address(N1));

    auto helloFromN2 = std::make_shared<inet::Olsrv2HelloPacket>();
    helloFromN2->setOriginator(inet::L3Address(N2));
    helloFromN2->setHTime(2.0);
    helloFromN2->setNeighAddrsArraySize(1);
    helloFromN2->setNeighAddrs(0, inet::L3Address(N1)); // 把 N1 放入邻居列表 => 对称判定
    nhdpN1.processHello(helloFromN2.get(), inet::L3Address(N2), now);

    // 从 N1 的 NHDP 数据生成 TC：应当通告其对称邻居 N2
    Olsrv2Core coreN1(N1);
    inet::Packet *tcPacket = coreN1.generateTc(nhdpN1.getDb(), now);
    expectTrue(tcPacket != nullptr, "generateTc must return a packet");

    auto tcChunk = tcPacket->peekAtFront<inet::Olsrv2TcGroup>();
    expectTrue(tcChunk != nullptr, "TC packet must contain Olsrv2TcGroup chunk");

    // 节点 S：拥有一跳对称邻居 N1，接收来自 N1 的 TC，计算到 N2 的路由
    Olsrv2Core coreS(S);
    NeighborTuple neighN1;
    neighN1.main_addr = N1;
    neighN1.is_symmetric = true;
    neighN1.validity.expires_at = now + (decoded_vtime_ms / 1000.0);
    coreS.state().upsertNeighbor(neighN1);

    coreS.processTcAndRecompute(tcChunk.get(), now);

    const auto& routes = coreS.state().getRoutes();
    const RouteTuple *rN2 = findRoute(routes, N2);
    expectTrue(rN2 != nullptr, "Route to N2 must exist after TC processing");
    if (rN2) {
        expectTrue(rN2->next_hop == N1, "Next hop to N2 must be N1");
        expectTrue(rN2->hop_count == 2, "Hop count to N2 should be 2 (S->N1->N2)");
    }

    delete tcPacket;
}

} // namespace

int run_integration_tests()
{
    std::cout << "Running test_integration...\n";
    testEndToEndHelloTcRoute();

    if (g_failures == 0) {
        std::cout << "[PASS] test_integration\n";
        return 0;
    }

    std::cerr << "[FAIL] test_integration failures=" << g_failures << "\n";
    return 1;
}

