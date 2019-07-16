import json

GeoJSON = json.load(open('AD-Map.geojson'))
items=GeoJSON['features']

R=[]
for i in range(len(items)):
	for j in range(len(items)):
		if 	i != j and \
			items[i]['properties']['p.dst'] == items[j]['properties']['p.org'] and	\
			items[i]['properties']['p.org'] == items[j]['properties']['p.dst']:
				if 	[items[i]['properties']['l.id'],items[j]['properties']['l.id']] not in R and \
					[items[j]['properties']['l.id'],items[i]['properties']['l.id']] not in R: 
					R.append([items[i]['properties']['l.id'],items[j]['properties']['l.id']])

print (len(R))
print (str(R).replace("], ","],\n"))