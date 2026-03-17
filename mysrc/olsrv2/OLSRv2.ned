//
// OLSRv2 routing protocol module (control plane over UDP).
//

package mysrc.olsrv2;

import inet.applications.contract.IApp;

@namespace(mysrc::olsrv2);
simple OLSRv2 like IApp
{
    parameters:
        @class(OLSRv2);
        @display("i=block/routing");

        string routingTableModule = default("^.ipv4.routingTable");
        string interfaceTableModule = default("^.interfaceTable");
        string networkProtocolModule = default("^.ipv4.ip");
        string mobilityModule = default("^.mobility");

        string interface = default("wlan0");
        string multicastGroup = default("224.0.0.109");
        int udpPort = default(698);

        double helloInterval @unit(s) = default(2s);
        double tcInterval @unit(s) = default(5s);
        double startJitter @unit(s) = default(0.2s);
    gates:
        input socketIn;
        output socketOut;
}
