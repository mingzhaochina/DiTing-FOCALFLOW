from obspy import read_inventory
import re
use_inv=True
if use_inv:
    invs = read_inventory("./stations/*.xml")
    with  open("stations.nll", "w") as o:
        for inv in invs:
            latitude = inv[0][0].latitude
            longitude = inv[0][0].longitude
            elevation =inv[0][0].elevation
            net = inv.code
            sta = inv[0].code
            channel=inv[0][0].code
            o.write(
                "{} {} {} {} {} {} {}\n".format("LOCSRCE",  sta, "LATLON", latitude,longitude, 0,elevation))
else:
    with  open("stations.nll", "w") as o,open ("../station.dat","r") as infile:
        f_csv = infile.readlines()
        for i in range(len(f_csv)):
            line = f_csv[i].strip('\n')
            lines = re.split('\s+', str(line))
            latitude = lines[1]
            longitude = lines[0]
            elevation = float(lines[5])*1000
            sta = lines[3]
            o.write(
                "{} {} {} {} {} {} {}\n".format("LOCSRCE", sta, "LATLON", latitude, longitude, 0, elevation))




