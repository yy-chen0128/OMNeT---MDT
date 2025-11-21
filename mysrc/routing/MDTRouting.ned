

package mysrc.routing;
@namespace(mysrc::routing);
simple MDTRouting{
    parameters:
    	//string routingTableModule = default("^.ipv4.routingTable");
        //string interfaceTableModule = default("^.interfaceTable");
        //string networkProtocolModule = default("^.ipv4.ip");
    	string routingTableModule;
        string interfaceTableModule;
        string networkProtocolModule;
        //string mdtDataModule;
        //double maintenanceInterval @unit(s);
        int udpPort;
        string interface;

        double activeRouteTimeout @unit(s);
    gates:
        input socketIn;
        output socketOut;

}