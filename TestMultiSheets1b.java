package tests;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;
import pack.M;
import pack.Pack;
import pack.Strip;
import processing.core.PApplet;

/**
 * 
 * @author Hao Hua, Southeast University Nanjing, http://labaaa.org
 *
 * 2D bin packing places a given set of polygons in standard single/multiple rectangular sheet(s), to minimize the use of the sheet(s).
 * This library does not involve other libraries, however, the example uses core.jar (https://processing.org) for a graphical interface.
 * The algorithm is effective when the ratio  (number of polygons / number of the types of polygons) is small.
 * 
 * 
 * 1. Input
 * a. Only simple polygon: no holes, no self-intersection.
 * b. Data structure: point - double[];  polygon - double[][]; all polygons - double[][][]. 
 * One might use other library to convert  other formats (e.g. dxf, obj) of polygon to double[][].
 * c. It's better to represent a polygon with a proper number of points. 
 * Too many points (e.g. a local detail contains dozens of points) slows down the algorithm. Avoid using too many points for a smooth curve.
 * One might use segment_max_length to create more points on a long edge of a polygon, if the polygon has very few points or the edge is very long.
 * d. The algorithm handles the "reference point" of a polygon internally, however, it is better to avoid the coordinates of points too far away from the origin. 
 * 
 * 
 * 2. Choose an algorithm
 * a. useAbey=true, Jostle heuristics for the 2D-irregular shapes bin packing problems with free rotation, R. P. Abeysooriya 2018
 * rotSteps:  Each polygon is rotated in the layout. Few steps (say 16) -> run fast  & poor result;  many steps (say 48) -> run slow & good results
 * segment_max_length: if (an edge of a polygon > segment_max_length),  breaks it into smaller segments. 
 * large value -> run fast  & poor result; small value -> run slow & good results
 * segment_max_length is related to translation steps.
 * 
 * b. useAbey=false, Waste minimization in irregular stock cutting, D. Dalalah, 2014
 * the rotation/translation steps depend on the polygons.
 * 
 * 
 * 3. Output
 * 
 * 
 * 
 * 
 */
public class TestMultiSheets1b extends PApplet {
	private DecimalFormat df=new DecimalFormat("##.###");
	private Random ran;
	
	//input
	private double[][][] randompolys;  //the input polygons
	private final float WID = 2400;  //width of the standard rectangular sheet
	private final float HEI = 1200;
	private final double margin= 6.0; 
	private final double preferX=0.499; // 0.501 or 1
	
	//parameters
	private final double segment_max_length =250.0; //250,400,800, use to break long edges if necessary, relative to the scale of the polgyons
	private final int rotSteps=36 ; //18,24,36,48, rotation steps
	private ArrayList<Pack> packs=new ArrayList<Pack>();
	private boolean useAbey=false;  
	// true: rotation steps depend on polygons,  R. P. Abeysooriya 2018
	// false: rotation steps are a prior,  D. Dalalah, 2014
	
	//output
	private int[] result_pack_id;  // result_pack_id[9]=2 means the 9th polygon is on the 2nd sheet.
	private double[][] result_cos_sin; // result_cos_sin[9]={0.5, 0.866} means the 9th polygon is rotated 60 degree w.r.t its reference point
	private double[][] result_position; // result_cos_sin[9] denotes the x,y-coordinate of the 9th polygon w.r.t its reference point

	public void setup() {
		size(1300, 800);
		int seed = 4346;
		ran = new Random(seed);
		randompolys = randomPolygons(210);// prepare the input polygons

		Double segment_len = useAbey ? null : segment_max_length;
		Pack pack = new Pack(randompolys, margin, segment_len, rotSteps, WID, HEI, preferX);
		pack.packOneSheet(useAbey);
		packs.add(pack);

		for (int i = 0; i < 100; i++) { // packing one sheet after another, 100 is estimated
			int size = packs.size();
			if (packs.get(size - 1).isEmpty()) {
				println(size + " sheets");
				break;
			}
			pack = packs.get(size - 1).clone();
			pack.packOneSheet(useAbey);
			packs.add(pack);
		}
		report();
	}

	private void report() {
		result_pack_id = new int[randompolys.length];
		result_cos_sin = new double[randompolys.length][];
		result_position = new double[randompolys.length][];
		for (int i = 0; i < packs.size(); i++) {
			Pack pack = packs.get(i);
			for (Strip strip : pack.fixs) {
				result_pack_id[strip.id] = i;
				result_cos_sin[strip.id] = strip.trigo;
				result_position[strip.id] = strip.position;
			}
		}
	}

	private double[][][] randomPolygons(int num) {
		double[][][] protos = new double[6][][];
		// protos[0] = new double[][] { { 0, 0 }, { 0, 2 }, { 6, 2 }, { 6, 12 }, { 8, 12 }, { 8, 0 } };
		protos[0] = new double[][] { { 0, 0 }, { 0, 10 }, { 10, 10 }, { 10, 0 } };
		protos[1] = new double[][] { { 0, 0 }, { 0, 16 }, { 1, 16 }, { 1, 0 } };
		protos[2] = new double[][] { { 0, 0 }, { 0, 12 }, { 9, 6 } };
		protos[3] = new double[][] { { 0, 0 }, { 0, 2 }, { 8, 2 }, { 8, 8 }, { 0, 8 }, { 0, 10 }, { 10, 10 }, { 10, 0 } };
		protos[4] = new double[][] { { 0, 0 }, { 0, 2 }, { 6, 2 }, { 6, 12 }, { 8, 12 }, { 8, 0 } };
		protos[5] = new double[][] { { 0, 0 }, { 0, 5 }, { 4, 1 }, { 8, 5 }, { 4, 10 }, { 0, 6 }, { 0, 11 }, { 9, 11 }, { 9, 0 } };

		double[][][] polys = new double[num][][];
		for (int i = 0; i < num; i++) {
			double[][] pl = protos[ran.nextInt(protos.length)];
			double sc = 8 + ran.nextDouble() * 12; // 8+ ran.nextDouble()*12 10+ ran.nextDouble()*20
			polys[i] = shake(3*sc, pl); //
		}
		for (int i = 0; i < num; i++) { 
			double[][] poly=polys[i] ;
			double dx=  ran.nextDouble() * 1000;
			double dy=  ran.nextDouble() * 1000;
			for(double[] p: poly){
				M._add(p, new double[]{dx,dy});
			}
		}
		return polys;
	}

	private double[][] shake(double s, double[][] ps) {
		double[][] arr = new double[ps.length][];
		double sng = ran.nextBoolean() ? s : (-s); // mirror
		double theta = ran.nextDouble() * 2 * PI;
		double cos = Math.cos(theta);
		double sin = Math.sin(theta);
		for (int i = 0; i < ps.length; i++) {
			double[] p = M.scale(sng, ps[i]);
			p[0] += (ran.nextDouble() - 0.5) * 0.25 * s;
			p[1] += (ran.nextDouble() - 0.5) * 0.25 * s;
			double x = cos * p[0] + sin * p[1];
			double y = -sin * p[0] + cos * p[1];
			arr[i] = new double[] { x, y };
		}
		return arr;
	}

	public void draw() {
		background(255);
		smooth();
		translate(40, 10);
		float sc = 0.15f;

		//display method 1
//		for (int i = 0; i < randompolys.length; i++) {
//			int pack_id = result_pack_id[i];
//			pushMatrix();
//			translate((pack_id % 3) * 380, (pack_id / 3) * 200); 
//			
//			noFill();
//			rect(0, 0, sc * WID, sc * HEI);
//			fill(0, 255, 0);
//			double[][] poly = randompolys[i];
//			poly = M.rotate(result_cos_sin[i], poly);
//			poly = M.move(result_position[i], poly);
//			draw(poly, sc);
//			
//			popMatrix();
//		}
		//display method 2
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 3; j++) {
				int id = i * 3 + j;
				if (id >= packs.size())
					return;
				Pack pack = packs.get(id);
				pushMatrix();
				translate(j * 380, i * 200); // (j * 270, i * 150
				noFill();
				rect(0, 0, sc * WID, sc * HEI);
				fill(0, 255, 0);
				for (Strip strip : pack.fixs) 
					draw(strip.inps, sc);
				popMatrix();
			}
		}
	}

	private void draw(double[][] ps, float sc) {
		beginShape();
		for (double[] p : ps)
			vertex((float) p[0] * sc, (float) p[1] * sc);
		endShape(CLOSE);
	}

}
