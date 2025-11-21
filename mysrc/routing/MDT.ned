
package mysrc.routing;

@namespace(mysrc::routing);
import inet.applications.contract.IApp;
import inet.networklayer.ipv4.Ipv4RoutingTable;
import inet.networklayer.common.InterfaceTable;

module MDT like IApp
{
    parameters:
        //MDT

        //@display("i=block/network2");
        @display("i=block/routing;is=s");

        double MaintenanceInterval @unit(s) = default(5s);
        double NSRPTimeout @unit(s) = default(5s);

        int coordinateDimension = default(3);
        double joinTimeout @unit(s) = default(5s);
		
		double activeRouteTimeout @unit(s) = default(20s); // the timeout value for cached routes If Hello messages are used


        int udpPort = default(654);
        string interface = default("wlan0");

        double keepAliveInterval @unit(s) = default(1s);
        double KANInterval @unit(s) = default(0.2s);
        double weaklinkKANInterval @unit(s) = default(0.15s);
        double neighborStaleTimeout @unit(s) = default(0.5s);
        double knownNodeStaleTimeout @unit(s) = default(15s);
        double tupleTimeout @unit(s) = default(10s);
        double dtInterval @unit(s) = default(2s);
        double initInterval @unit(s) = default(5s);

    gates:
        input socketIn;
        output socketOut;
    submodules:
        core: MDTRouting {
            parameters:
                @class("MDTRouting");
                udpPort = parent.udpPort;
                interface = parent.interface;
                activeRouteTimeout = parent.activeRouteTimeout;
                @display("i=block/control;x=0,y=0");
        }
        joinProtocol: JoinProtocol {
            parameters:
                @class("JoinProtocol");
                @display("i=block/protocol;x=100,y=50");
        }
        forwardProtocol: ForwardProtocol {
            parameters:
                @class("ForwardProtocol");
                @display("i=block/protocol;x=100,y=150");
        }
        maintenanceProtocol: MaintenanceProtocol {
            parameters:
                @class("MaintenanceProtocol");
                MaintenanceInterval = parent.MaintenanceInterval;
                NSRPTimeout = parent.NSRPTimeout;
                @display("i=block/protocol;x=150,y=250");
        }
		neighbor:Neighbor{
			parameters:
				@class("Neighbor");
				@display("i=block/protocol;x=100,y=0");
                keepAliveInterval = parent.keepAliveInterval;
                KANInterval = parent.KANInterval;
                weaklinkKANInterval = parent.weaklinkKANInterval;
                neighborStaleTimeout = parent.neighborStaleTimeout;
                knownNodeStaleTimeout = parent.knownNodeStaleTimeout;
                tupleTimeout = parent.tupleTimeout;
                dtInterval = parent.dtInterval;
                initInterval = parent.initInterval;
		}
	connections:
	    socketIn --> core.socketIn;
	    core.socketOut --> socketOut;
}

