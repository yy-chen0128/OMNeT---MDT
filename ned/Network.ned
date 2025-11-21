
package ned;

import inet.networklayer.configurator.ipv4.Ipv4NetworkConfigurator;
import inet.physicallayer.wireless.ieee80211.packetlevel.Ieee80211ScalarRadioMedium;
import inet.visualizer.common.IntegratedVisualizer;
import inet.environment.common.PhysicalEnvironment;
import inet.node.aodv.AodvRouter;
import mysrc.RPC.OmnetRpcLogger;

network DynamicAdHocNetwork
{
    parameters:
        int numNodes = default(25);
        double playgroundSizeX @unit(m) = default(1500m);
        double playgroundSizeY @unit(m) = default(1500m);
        double playgroundSizeZ @unit(m) = default(1500m);
        @display("bgb=1500,1500");
    submodules:
        //coordinateSystem: GlobalCoordinateSystem {
          //  @display("p=50,250");
        //}
        visualizer: IntegratedVisualizer{
        	parameters:
                @display("p=50,50");
        }
        configurator: Ipv4NetworkConfigurator {
            parameters:
                @display("p=101,50");
                // allocate address
                addressing.typename = "FlatNetworkConfigurator";
                config = xml("<config> \
                    <interface hosts='**' address='10.0.0.x' netmask='255.255.255.0'/> \
                </config>");
        }
		physicalEnvironment: PhysicalEnvironment {
            parameters:
                @display("p=300,50");
        }
        radioMedium: Ieee80211ScalarRadioMedium {
            parameters:
                @display("p=50,100");
        }

		rpcLogger: OmnetRpcLogger {
            parameters:
                @display("p=20,20");
                rpcPort = default(9002);
                rpcPollInterval = default(0.1s);
        }

        node[numNodes]: Node {
            parameters:
                //@display("p= , m,50;is=s");
        }
}



