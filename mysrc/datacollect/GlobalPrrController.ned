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

package mysrc.datacollect;

simple GlobalPrrController
{
    parameters:
        string mode = default("fast");    // "fast" or "physical"
        @display("i=block/cog");

		// SIGNALS
        @signal[reachable](type=bool);
        @signal[receiverCount](type=int);
        @signal[theoreticalPRR](type=double);

        // STATISTICS (note: use semicolons inside the parentheses and terminate each @statistic with semicolon)
        @statistic[reachable](title="Theoretical reachability"; source=reachable; record=vector);
        @statistic[receiverCount](title="Theoretical receiver count"; source=receiverCount; record=vector);
        @statistic[theoreticalPRR](title="Theoretical PRR per packet"; source=theoreticalPRR; record=vector);
}
