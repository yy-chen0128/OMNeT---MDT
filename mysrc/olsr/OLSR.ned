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

package mysrc.olsr;
import inet.applications.contract.IApp;
import inet.networklayer.ipv4.Ipv4RoutingTable;
import inet.networklayer.common.InterfaceTable;
@namespace(mysrc::olsr);
simple OLSR like IApp
{
	parameters:
        @class(OLSR);
        @display("i=block/network2");
        string routingTableModule;
        string interfaceTableModule;
        string networkProtocolModule;
        string interface = default("wlan0");

        // OLSR specific parameters
        //int Willingness @unit(s) = default(3s);  // node's willingness to carry traffic for others
        double helloInterval @unit(s) = default(0.75s);  // interval for HELLO messages in seconds
        double tcInterval @unit(s) = default(1.5s);     // interval for Topology Control (TC) messages in seconds
        double tcCacheInterval @unit(s) = default(0.05s);

    gates:
        input socketIn;
        output socketOut;
}
