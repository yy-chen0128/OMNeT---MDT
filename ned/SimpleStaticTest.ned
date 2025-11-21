//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see http://www.gnu.org/licenses/.
//

package ned;

import inet.networklayer.configurator.ipv4.Ipv4NetworkConfigurator;
import inet.physicallayer.wireless.ieee80211.packetlevel.Ieee80211ScalarRadioMedium;
import inet.visualizer.common.IntegratedVisualizer;
import ned.Node;

network SimpleStaticTest
{
    parameters:
        @display("bgb=800,800");
    submodules:
        visualizer: IntegratedVisualizer {
            parameters: @display("p=50,50");
        }
        configurator: Ipv4NetworkConfigurator {
            parameters:
                @display("p=120,50");
                addressing.typename = "FlatNetworkConfigurator";
                config = xml("<config><interface hosts='**' address='10.0.0.x' netmask='255.255.255.0'/></config>");
        }
        radioMedium: Ieee80211ScalarRadioMedium {
            parameters: @display("p=50,150");
        }

        node[7]: Node {
            parameters:
                //@display("p=200,200");
        }
}



