
named_points = {
    'SCCI' : [-53.01062,-70.85168,42],
    'SCAR' : [-18.3483,-70.3386, 167],
    'SCTE' : [-41.438611, -73.0939, 294],
    'SCVD' : [-39.649722,-73.086111,59]}

airport_lon = []
airport_lat = []
for k,v in named_points.items():
    airport_lon.append(v[1])
    airport_lat.append(v[0])
