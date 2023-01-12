import fiona
import shapely
from shapely.ops import cascaded_union

sellist =  ['Región de Valparaíso',
            'Región Metropolitana de Santiago',
            'Región del Maule',
            "Región del Libertador Bernardo O'Higgins"]

fname = '../../../megafires_data/shp/Regiones/Regional.shp'

shapes = list()
with fiona.open(fname) as source:
    for feature in source:
        geom = feature['geometry']
        prop = feature['properties']
        # print(len(geom['coordinates']), prop['Region'])
        if prop['Region'] in sellist:
            shp = shapely.geometry.shape(geom)
            if prop['Region'] == 'Región de Valparaíso':
                shapes.extend([x for x in shp.geoms if x.is_valid])
            elif prop['Region'] == "Región del Libertador Bernardo O'Higgins":
                shapes.append(shp)
            else:
                shapes.append(shp)
shapes = cascaded_union(shapes)


