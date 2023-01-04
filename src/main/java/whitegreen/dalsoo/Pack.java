package whitegreen.dalsoo;

import static java.lang.Math.PI;

import java.util.ArrayList;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import net.jafama.FastMath;

/**
 * 
 * @author Hao Hua, Southeast University, whitegreen@163.com
 * @author Michael Carleton
 *
 */
public final class Pack {

	private static final GeometryFactory GEOM_FACT = new GeometryFactory();

	public ArrayList<Strip> fixs = new ArrayList<>();
	public ArrayList<Strip> movs = new ArrayList<>();
	public Convex cntConvex;
	/**
	 * Sin and cos values for each rotation step.
	 */
	private final double[][] trigos;
	private final int rotSteps;
	private final double width, height;
	private static final double areasc = 1E-6; //
	private final double preferX; // 0.501 or 1

	Pack(double[][] trigos, int rotSteps, double width, double height, double preferX) {
		this.trigos = trigos;
		this.rotSteps = rotSteps;
		this.width = width;
		this.height = height;
		this.preferX = preferX;
	}

	/**
	 * 
	 * @param polys              array of polygons [p1, p2]; for each polygon, its
	 *                           vertices are expressed in [x, y] coordinate pairs:
	 *                           [[p1x, p1y], [p2x, p2y]...]
	 * @param spacing            boundary spacing between packed objects
	 * @param segment_max_length
	 * @param rotSteps
	 * @param width              width of each bin
	 * @param height             height of each bin
	 * @param preferX
	 */
	public Pack(double[][][] polys, double spacing, Double segment_max_length, int rotSteps, double width, double height,
			double preferX) {
		if (0 > preferX || 1 < preferX) {
			throw new IllegalArgumentException("skew must be between 0 and 1 (inclusive)");
		}
		for (int i = 0; i < polys.length; i++) {
			double[][] poly = polys[i];
			Strip strip = new Strip(i, poly, spacing, segment_max_length);
			movs.add(strip);
		}
//		Collections.sort(movs); // sort by area, smallest first
		this.preferX = preferX;

		this.rotSteps = rotSteps;
		trigos = new double[rotSteps][];
		for (int i = 0; i < rotSteps; i++) {
			double theta = i * 2 * PI / rotSteps;
			trigos[i] = new double[] { FastMath.cos(theta), FastMath.sin(theta) };
		}
		this.width = width;
		this.height = height;
	}

	public void packOneSheet(boolean abey) {
		place1stStrip();
		int size = movs.size();
		for (int i = 0; i < size; i++) {
			// place next largest polygon
			if (abey) {
				placeAnotherStrip_Abey(size - 1 - i);
			} else {
				placeAnotherStrip_Dalalah(size - 1 - i);
			}
		}
		ArrayList<Strip> list = new ArrayList<>();
		for (Strip stp : movs) {
			if (null == stp) {
				continue;
			}
			list.add(stp);
		}
		movs = list;
	}

	private void place1stStrip() {
		int rotid = 0;
		double minArea = 1000000000;
		int sid = movs.size() - 1; // ***************************************** last one
		Strip first = movs.get(sid);
		for (int i = 0; i < rotSteps; i++) {
			double[][] tp = M.rotate(trigos[i], first.outps);
			double[] bd = M.boundBox(tp); // minx, maxx, miny, maxy
			if (bd[1] - bd[0] > width || bd[3] - bd[2] > height) {
				continue;
			}
			double area = areasc * (bd[1] - bd[0]) * (bd[3] - bd[2]);
			double[] center = M.mean(tp);
			double len = preferX * (center[0] - bd[0]) + (1 - preferX) * (center[1] - bd[2]);
			area *= len;
			if (minArea > area) {
				minArea = area;
				rotid = i;
			}
		}
		double[][] tp = M.rotate(trigos[rotid], first.outps);
		double[] bd = M.boundBox(tp);
		first.fix_rotate_move(trigos[rotid], new double[] { -bd[0], -bd[2] });
		movs.remove(sid);
		placeStrip(first);
		cntConvex = new Convex(first);
	}

	private boolean placeAnotherStrip_Abey(int sid) {
		Strip stp = movs.get(sid);
		double min = 1000000000;
		double[] min_cossin = null;
		double[] min_trans = null;
		Convex min_con = null;
		double[][] opl = stp.outps;
		for (int i = 0; i < opl.length; i++) { // each vertex of new strip
			double[] p = opl[i];
			double[] d0 = M.sub(opl[(i - 1 + opl.length) % opl.length], p);
			double[] d2 = M.sub(opl[(i + 1) % opl.length], p);
			double mag0 = M.mag(d0);
			double mag2 = M.mag(d2);
			double cb0 = d0[0] / mag0;
			double sb0 = d0[1] / mag0;
			double cb2 = d2[0] / mag2;
			double sb2 = d2[1] / mag2;

			for (Strip fixed : fixs) {
				double[][] fopl = fixed.outps;
				for (int j = 0; j < fopl.length; j++) { // each vertex of each fixed strip
					double[] v = fopl[j];
					double[] t0 = M.sub(fopl[(j - 1 + fopl.length) % fopl.length], v);
					double[] t2 = M.sub(fopl[(j + 1) % fopl.length], v);
					double m0 = M.mag(t0);
					double m2 = M.mag(t2);
					double ca0 = t0[0] / m0;
					double sa0 = t0[1] / m0;
					double ca2 = t2[0] / m2;
					double sa2 = t2[1] / m2;

					for (int h = 0; h < 2; h++) { // two angles
						double[] cossin;
						if (0 == h) {
							cossin = new double[] { ca0 * cb2 + sa0 * sb2, sa0 * cb2 - ca0 * sb2 }; // a0 - b2
						} else {
							cossin = new double[] { ca2 * cb0 + sa2 * sb0, sa2 * cb0 - ca2 * sb0 };// a2 - b0
						}

						double[][] rot_opl = M.rotate(cossin, opl);
						double[] trans = M.sub(v, rot_opl[i]); // *****************
						double[][] trans_rot_outpoly = M.move(trans, rot_opl);
						if (feasible(trans_rot_outpoly)) {
							double[][] trans_rot_inpoly = M.move(trans, M.rotate(cossin, stp.inps));
							Convex tmpcon = cntConvex.clone();
							for (double[] trp : trans_rot_inpoly) {
								tmpcon.increment_hull(trp);
							}

							double conarea = areasc * M.areaAbs(tmpcon.convex);
							double[] center = M.mean(trans_rot_inpoly);
							double area;
							area = conarea * (preferX * Math.abs(center[0]) + (1 - preferX) * Math.abs(center[1]));
							if (min > area) {
								min = area;
								min_cossin = cossin;
								min_trans = trans;
								min_con = tmpcon;
							}
						}
					} // for h
				}
			}
		} // for each vertex of new strip
		if (null == min_cossin) { // no solution, stop
			return false;
		}
		stp.fix_rotate_move(min_cossin, min_trans);
		movs.set(sid, null);
		placeStrip(stp);
		cntConvex = min_con;
		return true;
	}

	private boolean placeAnotherStrip_Dalalah(int sid) {
		Strip stp = movs.get(sid);
		double min = 1000000000;
		int min_rotid = -1;
		double[] min_trans = null;
		Convex min_con = null;
		for (int i = 0; i < rotSteps; i++) {// each angle of new strip
			double[][] rotated_outpoly = M.rotate(trigos[i], stp.outps);
			double[][] rotated_inpoly = M.rotate(trigos[i], stp.inps);

			for (double[] p : rotated_outpoly) { // each vertex of new strip
				for (Strip fixed : fixs) {
					for (double[] v : fixed.outps) { // each vertex of each fixed strip
						double[] trans = M.sub(v, p);
						double[][] trans_rot_outpoly = M.move(trans, rotated_outpoly);
						if (feasible(trans_rot_outpoly)) {
							double[][] trans_rot_inpoly = M.move(trans, rotated_inpoly);
							Convex tmpcon = cntConvex.clone();
							for (double[] trp : trans_rot_inpoly) {
								tmpcon.increment_hull(trp);
							}

							double conarea = areasc * M.areaAbs(tmpcon.convex);
							double[] center = M.mean(trans_rot_outpoly);// trans_rot_outpoly
							double area = conarea * (preferX * Math.abs(center[0]) + Math.abs(center[1]));
							if (min > area) {
								min = area;
								min_rotid = i;
								min_trans = trans;
								min_con = tmpcon;
							}
						}
					}
				}
			} // for each vertex of new strip
		} // for each angle of new strip
		if (0 > min_rotid) { // no solution, stop
			return false;
		}
		stp.fix_rotate_move(trigos[min_rotid], min_trans);
		// movs.remove(sid);
		movs.set(sid, null);
		placeStrip(stp);
		cntConvex = min_con;
		return true;
	}

	private boolean feasible(double[][] poly) {
		for (double[] p : poly) {
			if (p[0] < 0 || p[0] > width || p[1] < 0 || p[1] > height) {
				// out of bounds -- infeasible!
				return false;
			}
		}
		final double[] bb = M.boundBox(poly);

		for (Strip fixed : fixs) {
			if (overlapFast(poly, bb, fixed)) {
				return false;
			}
		}
		return true;
	}

	/**
	 * @deprecated
	 */
	private static final boolean overlap(final double[][] poly1, final double[][] poly2) {
		if (!M.intersect_boundBox(poly1, poly2)) {
			return false;
		}
		double[][] sml, big;
		if (M.areaAbs(poly1) < M.areaAbs(poly2)) {
			sml = poly1;
			big = poly2;
		} else {
			sml = poly2;
			big = poly1;
		}
		for (double[] p : sml) {
			if (M.inside(p, big)) {
				return true;
			}
		}
		for (int i = 0; i < sml.length; i++) {
			for (int j = 0; j < big.length; j++) {
				if (M.intersection(sml[i], sml[(i + 1) % sml.length], big[j], big[(j + 1) % big.length])) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Determines whether poly overlaps with strip.
	 * 
	 * @param poly1     vertices of poly 1
	 * @param poly1BB   bounding box of poly 1
	 * @param poly1Area area of poly 1
	 * @param poly2     vertices of poly 2
	 * @return
	 */
	private static final boolean overlapFast(final double[][] poly1, final double[] poly1BB, final Strip strip) {

		// 1. test bounding boxes
		if (!M.intersect_boundBoxes(poly1BB, strip.bb)) {
			return false;
		}
		// 2. test the centroid against candidate poly
		/*
		 * NOTE this isn't completely robust (a robust method would test every point of
		 * strip against candidate) but would catch almost all occurrences (so is worth
		 * the speed-up overall).
		 */
		if (M.inside(strip.centroid, poly1)) {
			return true;
		}

		// 3. test edge intersections
		for (int i = 0; i < poly1.length; i++) {
			for (int j = 0; j < strip.inps.length; j++) {
				if (M.intersectionFast(poly1[i], poly1[(i + 1) % poly1.length], strip.inps[j],
						strip.inps[(j + 1) % strip.inps.length])) {
					return true;
				}
			}
		}
		return false;

	}

	private void placeStrip(Strip strip) {
		strip.fix();
		fixs.add(strip);
	}

	static Polygon toPolygon(double[][] poly) {
		Coordinate[] coords = new Coordinate[poly.length + 1];
		for (int i = 0; i < coords.length - 1; i++) {
			coords[i] = new Coordinate(poly[i][0], poly[i][1]);
		}
		coords[coords.length - 1] = new Coordinate(poly[0][0], poly[0][1]); // close
		return GEOM_FACT.createPolygon(coords);
	}

	static double[][] toArray(Coordinate[] coords) {
		double[][] poly = new double[coords.length - 1][2];
		for (int i = 0; i < coords.length - 1; i++) { // -1, unclose
			poly[i][0] = coords[i].x;
			poly[i][1] = coords[i].y;
		}
		return poly;
	}

	public double leftOverArea() {
		double sum = 0;
		for (Strip strip : movs) {
			sum += strip.inarea;
		}
		return sum;
	}

	public int[] getStripIds() {
		int size = fixs.size();
		int[] ids = new int[size];
		for (int i = 0; i < size; i++) {
			ids[i] = fixs.get(i).id;
		}
		return ids;
	}

	public double[][] getStripRotations() {
		int size = fixs.size();
		double[][] rots = new double[size][];
		for (int i = 0; i < size; i++) {
			rots[i] = fixs.get(i).trigo;
		}
		return rots;
	}

	public double[][] getStripPositions() {
		int size = fixs.size();
		double[][] ps = new double[size][];
		for (int i = 0; i < size; i++) {
			ps[i] = fixs.get(i).position;
		}
		return ps;
	}

	@Override
	public Pack clone() {
		Pack np = new Pack(trigos, rotSteps, width, height, preferX);
		np.movs = movs;
		np.cntConvex = cntConvex;
		return np;
	}

	public boolean isEmpty() {
		return movs.isEmpty();
	}

}
