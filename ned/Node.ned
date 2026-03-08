package ned;

import inet.node.inet.AdhocHost;
import mysrc.routing.MDT;

module Node extends AdhocHost
{
    parameters:
        @display("i=block/routing");
        //*.networkLayer.routing.typename = "MDTRouting";
        *.routing.typename = "MDT";
        //*.routingApp.typename = "MDTRouting";
     submodules:
        MDTrouting: MDT {
           parameters:
                @display("p=825,226");
       }
    connections:
        MDTrouting.socketOut --> at.in++;
        MDTrouting.socketIn <-- at.out++;
}



