# Dalsoo-Bin-Packing
Bin packing places a given set of polygons in standard single/multiple rectangular sheet(s), to minimize the use of the sheet(s).
 This library does not involve other libraries, however, the example uses core.jar (https://processing.org) for a graphical interface.
 The algorithm is effective when the ratio  (number of polygons / number of the types of polygons) is small.

![alt text](relative/path/to/img.jpg?raw=true "multiple sheets")
a
![alt text](/multiple sheets.jpg "Description goes here")

1. Input 

a. Only simple polygon: no holes, no self-intersection.

b. Data structure: point - double[];  polygon - double[][]; all polygons - double[][][]. 
One might use other library to convert  other formats (e.g. dxf, obj) of polygon to double[][].

c. It's better to represent a polygon with a proper number of points. 
Too many points (e.g. a local detail contains dozens of points) slows down the algorithm. Avoid using too many points for a smooth curve.
One might use segment_max_length to create more points on a long edge of a polygon, if the polygon has very few points or the edge is very long.

d. The algorithm handles the "reference point" of a polygon internally, however, it is better to avoid the coordinates of points too far away from the origin. 



2. Choose an algorithm

a. useAbey=true, Jostle heuristics for the 2D-irregular shapes bin packing problems with free rotation, R. P. Abeysooriya 2018
rotSteps:  Each polygon is rotated in the layout. Few steps (say 16) -> run fast  & poor result;  many steps (say 48) -> run slow & good results
segment_max_length: if (an edge of a polygon > segment_max_length),  breaks it into smaller segments. 
large value -> run fast  & poor result; small value -> run slow & good results
segment_max_length is related to translation steps.

b. useAbey=false, Waste minimization in irregular stock cutting, D. Dalalah, 2014
the rotation/translation steps depend on the polygons.



3. Output options

a. obtain the transformation of each polygon in the final layout, including:
int[] result_pack_id;  // result_pack_id[9]=2 means the 9th polygon is on the 2nd sheet.
double[][] result_cos_sin; // result_cos_sin[9]={0.5, 0.866} means the 9th polygon is rotated 60 degree w.r.t its reference point
double[][] result_position; // result_cos_sin[9] denotes the x,y-coordinate of the 9th polygon w.r.t its reference point

b. obtain the geometry of each in the final layout 
Pack pack = packs.get(id);   //obtain a sheet (pack) by id.
for (Strip strip : pack.fixs) { // obtain a placed strip (polygon) in a sheet.
      strip.inps     // the polgyon's points, one might use other library to convert double[][] to other formats (dxf, obj) of polygon 
}
