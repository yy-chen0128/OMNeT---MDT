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

package mysrc.routing;
@namespace(mysrc::routing);
simple Neighbor{
    parameters:
        @display("i=block/protocol");
        double keepAliveInterval @unit(s);
        double KANInterval @unit(s);
        double weaklinkKANInterval @unit(s);
        double neighborStaleTimeout @unit(s);
        double knownNodeStaleTimeout @unit(s);
        double tupleTimeout @unit(s);
        double dtInterval @unit(s);
        double initInterval @unit(s);
}

