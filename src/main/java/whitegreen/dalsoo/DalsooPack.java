package whitegreen.dalsoo;

import static java.lang.Math.PI;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

import net.jafama.FastMath;

/**
 * Bin-packs irregular shapes sequentially.
 * <p>
 * Algorithms are adapted from the papers of Abeysooriya 2018 and Dalalah 2014.
 * Bin packing places a given set of polygons in standard single/multiple
 * rectangular sheet(s), minimising the use of the sheet(s).
 * 
 * @author Hao Hua, Southeast University, whitegreen@163.com
 * @author Michael Carleton
 *
 */
public final class DalsooPack {

	private static final double AREA_SC = 1E-6;

	public List<PackedPoly> packedPolys = new ArrayList<>();
	List<PackedPoly> pendingPolys = new ArrayList<>();
	Convex cntConvex;

	/** Sin and cos values for each rotation step. */
	private final double[][] trigos;
	private final int rotSteps;
	/** bin dimensions */
	private final double binWidth, binHeight;
	private final double preferX;

	DalsooPack(double[][] trigos, int rotSteps, double width, double height, double preferX) {
		this.trigos = trigos;
		this.rotSteps = rotSteps;
		this.binWidth = width;
		this.binHeight = height;
		this.preferX = preferX;
	}

	/**
	 * Creates an packing. Call {@link #packOneBin(boolean, boolean) packOneBin()}
	 * to start packing.
	 * 
	 * @param polys              array of polygons [p1, p2]; for each polygon, its
	 *                           vertices are expressed in [x, y] coordinate pairs:
	 *                           [[p1x, p1y], [p2x, p2y]...]. Polygons should be
	 *                           simple, having no holes and no self-intersection.
	 * @param spacing            boundary spacing between packed objects.
	 * @param segment_max_length
	 * @param rotSteps           number of unique rotation variants of each polygon.
	 *                           a larger number results in higher-quality packing
	 *                           but longer runtime (not used in abey packing).
	 * @param width              width of each bin
	 * @param height             height of each bin
	 * @param hSkew              horizontal skew [0...1]. Determines to what extent
	 *                           shapes are packed horizontally (vs vertically).
	 *                           i.e. When =1.0, shapes will start packing against
	 *                           the y axis, and fill left-to-right (horizontal);
	 *                           When =0.0, shapes will pack from the x-axis and
	 *                           fill top-to-bottom (vertical).
	 */
	public DalsooPack(double[][][] polys, double spacing, Double segment_max_length, int rotSteps, double width, double height,
			double hSkew) {
		if (0 > hSkew || 1 < hSkew) {
			throw new IllegalArgumentException("skew must be between 0 and 1 (inclusive)");
		}
		for (int i = 0; i < polys.length; i++) {
			double[][] poly = polys[i];
			PackedPoly p = new PackedPoly(i, poly, spacing, segment_max_length);
			pendingPolys.add(p);
		}
		this.preferX = hSkew;

		this.rotSteps = rotSteps;
		trigos = new double[rotSteps][];
		for (int i = 0; i < rotSteps; i++) {
			double theta = i * 2 * PI / rotSteps;
			trigos[i] = new double[] { FastMath.cos(theta), FastMath.sin(theta) };
		}
		this.binWidth = width;
		this.binHeight = height;
	}

	/**
	 * Packs all polygons, possibly across multiple bins.
	 * 
	 * @param abey whether to use
	 * @return
	 */
	public List<DalsooPack> packAll(boolean abey) {
		/*
		 * The pack essentially forks itself with the remaining polygons to model new
		 * bins.
		 */
		List<DalsooPack> packedBins = new ArrayList<>();
		DalsooPack pack = this;
		do {
			pack.packOneBin(abey, false);
			packedBins.add(pack);
			pack = pack.clone();
		} while (!packedBins.get(packedBins.size() - 1).isEmpty());

		return packedBins;
	}

	/**
	 * 
	 * @param abey
	 * @param largestFirst whether to sort the polygons by area (largest first)
	 *                     before packing
	 */
	public void packOneBin(boolean abey, boolean largestFirst) {
		if (largestFirst) {
			Collections.sort(pendingPolys); // sort by area, smallest first
		}
		packFirstPoly();
		int size = pendingPolys.size();
		for (int i = 0; i < size; i++) {
			// place next largest polygon
			if (abey) {
				packPolyAbey(size - 1 - i);
			} else {
				packPolyDalalah(size - 1 - i);
			}
		}

		pendingPolys.removeIf(Objects::isNull);
	}

	private void packFirstPoly() {
		int rotid = 0;
		double minArea = Double.MAX_VALUE;
		int sid = pendingPolys.size() - 1; // ***************************************** last one
		PackedPoly first = pendingPolys.get(sid);
		for (int i = 0; i < rotSteps; i++) {
			double[][] tp = MathUtil.rotate(trigos[i], first.outpts);
			double[] bd = MathUtil.boundBox(tp); // minx, maxx, miny, maxy
			if (bd[1] - bd[0] > binWidth || bd[3] - bd[2] > binHeight) {
				continue;
			}
			double area = AREA_SC * (bd[1] - bd[0]) * (bd[3] - bd[2]);
			double[] center = MathUtil.mean(tp);
			double len = preferX * (center[0] - bd[0]) + (1 - preferX) * (center[1] - bd[2]);
			area *= len;
			if (minArea > area) {
				minArea = area;
				rotid = i;
			}
		}
		double[][] tp = MathUtil.rotate(trigos[rotid], first.outpts);
		double[] bd = MathUtil.boundBox(tp);
		first.fix_rotate_move(trigos[rotid], new double[] { -bd[0], -bd[2] });
		pendingPolys.remove(sid);
		placePackedPoly(first);
		cntConvex = new Convex(first);
	}

	private boolean packPolyAbey(int sid) {
		PackedPoly poly = pendingPolys.get(sid);
		double minArea = Double.MAX_VALUE;
		double[] min_cossin = null;
		double[] min_trans = null;
		Convex min_con = null;
		double[][] opl = poly.outpts;
		for (int i = 0; i < opl.length; i++) { // each vertex of new poly
			double[] p = opl[i];
			double[] d0 = MathUtil.sub(opl[(i - 1 + opl.length) % opl.length], p);
			double[] d2 = MathUtil.sub(opl[(i + 1) % opl.length], p);
			double mag0 = MathUtil.mag(d0);
			double mag2 = MathUtil.mag(d2);
			double cb0 = d0[0] / mag0;
			double sb0 = d0[1] / mag0;
			double cb2 = d2[0] / mag2;
			double sb2 = d2[1] / mag2;

			for (PackedPoly fixed : packedPolys) {
				double[][] fopl = fixed.outpts;
				for (int j = 0; j < fopl.length; j++) { // each vertex of each fixed poly
					double[] v = fopl[j];
					double[] t0 = MathUtil.sub(fopl[(j - 1 + fopl.length) % fopl.length], v);
					double[] t2 = MathUtil.sub(fopl[(j + 1) % fopl.length], v);
					double m0 = MathUtil.mag(t0);
					double m2 = MathUtil.mag(t2);
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

						double[][] rot_opl = MathUtil.rotate(cossin, opl);
						double[] trans = MathUtil.sub(v, rot_opl[i]); // *****************
						double[][] trans_rot_outpoly = MathUtil.move(trans, rot_opl);
						if (isFeasible(trans_rot_outpoly)) {
							double[][] trans_rot_inpoly = MathUtil.move(trans, MathUtil.rotate(cossin, poly.inpts));
							Convex tmpcon = cntConvex.clone();
							for (double[] trp : trans_rot_inpoly) {
								tmpcon.increment_hull(trp);
							}

							double conarea = AREA_SC * MathUtil.areaAbs(tmpcon.convex);
							double[] center = MathUtil.mean(trans_rot_inpoly);
							double area;
							area = conarea * (preferX * Math.abs(center[0]) + (1 - preferX) * Math.abs(center[1]));
							if (minArea > area) {
								minArea = area;
								min_cossin = cossin;
								min_trans = trans;
								min_con = tmpcon;
							}
						}
					} // for h
				}
			}
		} // for each vertex of new poly
		if (null == min_cossin) { // no solution, stop
			return false;
		}
		poly.fix_rotate_move(min_cossin, min_trans);
		pendingPolys.set(sid, null); // done, mark as null
		placePackedPoly(poly);
		cntConvex = min_con;
		return true;
	}

	private boolean packPolyDalalah(int sid) {
		PackedPoly stp = pendingPolys.get(sid);
		double minArea = Double.MAX_VALUE;
		int min_rotid = -1;
		double[] min_trans = null;
		Convex min_con = null;
		for (int i = 0; i < rotSteps; i++) {// each angle of new poly
			double[][] rotated_outpoly = MathUtil.rotate(trigos[i], stp.outpts);
			double[][] rotated_inpoly = MathUtil.rotate(trigos[i], stp.inpts);

			for (double[] p : rotated_outpoly) { // each vertex of new poly
				for (PackedPoly fixed : packedPolys) {
					for (double[] v : fixed.outpts) { // each vertex of each fixed poly
						double[] trans = MathUtil.sub(v, p);
						double[][] trans_rot_outpoly = MathUtil.move(trans, rotated_outpoly);
						if (isFeasible(trans_rot_outpoly)) {
							double[][] trans_rot_inpoly = MathUtil.move(trans, rotated_inpoly);
							Convex tmpcon = cntConvex.clone();
							for (double[] trp : trans_rot_inpoly) {
								tmpcon.increment_hull(trp);
							}

							double conarea = AREA_SC * MathUtil.areaAbs(tmpcon.convex);
							double[] center = MathUtil.mean(trans_rot_outpoly);// trans_rot_outpoly
							double area = conarea * (preferX * Math.abs(center[0]) + Math.abs(center[1]));
							if (minArea > area) {
								minArea = area;
								min_rotid = i;
								min_trans = trans;
								min_con = tmpcon;
							}
						}
					}
				}
			} // for each vertex of new poly
		} // for each angle of new poly
		if (0 > min_rotid) { // no solution, stop
			return false;
		}
		stp.fix_rotate_move(trigos[min_rotid], min_trans);
		// movs.remove(sid);
		pendingPolys.set(sid, null);
		placePackedPoly(stp);
		cntConvex = min_con;
		return true;
	}

	/**
	 * Determines whether the potential polygon can be placed.
	 */
	private boolean isFeasible(double[][] poly) {
		for (double[] p : poly) {
			if (p[0] < 0 || p[0] > binWidth || p[1] < 0 || p[1] > binHeight) {
				// out of bounds -- infeasible!
				return false;
			}
		}
		final double[] bb = MathUtil.boundBox(poly);

		for (PackedPoly fixed : packedPolys) {
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
		if (!MathUtil.intersect_boundBox(poly1, poly2)) {
			return false;
		}
		double[][] sml, big;
		if (MathUtil.areaAbs(poly1) < MathUtil.areaAbs(poly2)) {
			sml = poly1;
			big = poly2;
		} else {
			sml = poly2;
			big = poly1;
		}
		for (double[] p : sml) {
			if (MathUtil.inside(p, big)) {
				return true;
			}
		}
		for (int i = 0; i < sml.length; i++) {
			for (int j = 0; j < big.length; j++) {
				if (MathUtil.intersection(sml[i], sml[(i + 1) % sml.length], big[j], big[(j + 1) % big.length])) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Determines whether potential poly overlaps with placed poly.
	 * 
	 * @param candidate   vertices of poly 1
	 * @param candidateBB bounding box of poly 1
	 * @param poly1Area   area of poly 1
	 * @param poly2       vertices of poly 2
	 * @return
	 */
	private static final boolean overlapFast(final double[][] candidate, final double[] candidateBB,
			final PackedPoly placedPoly) {

		// 1. test bounding boxes
		if (!MathUtil.intersect_boundBoxes(candidateBB, placedPoly.bb)) {
			return false;
		}
		// 2. test the centroid against candidate poly
		/*
		 * NOTE this isn't completely robust (a robust method would test every point of
		 * placed poly against candidate) but will catch almost all occurrences (so is
		 * worth the speed-up overall).
		 */
		if (MathUtil.inside(placedPoly.centroid, candidate)) {
			return true;
		}

		// 3. test edge intersections
		for (int i = 0; i < candidate.length; i++) {
			for (int j = 0; j < placedPoly.inpts.length; j++) {
				if (MathUtil.intersectionFast(candidate[i], candidate[(i + 1) % candidate.length], placedPoly.inpts[j],
						placedPoly.inpts[(j + 1) % placedPoly.inpts.length])) {
					return true;
				}
			}
		}
		return false;

	}

	private void placePackedPoly(PackedPoly poly) {
		poly.place();
		packedPolys.add(poly);
	}

	public double getEmptyArea() {
		double sum = 0;
		for (PackedPoly poly : pendingPolys) {
			sum += poly.inarea;
		}
		return sum;
	}

	public int[] getPackingIds() {
		int size = packedPolys.size();
		int[] ids = new int[size];
		for (int i = 0; i < size; i++) {
			ids[i] = packedPolys.get(i).id;
		}
		return ids;
	}

	public double[][] getPackingRotations() {
		int size = packedPolys.size();
		double[][] rots = new double[size][];
		for (int i = 0; i < size; i++) {
			rots[i] = packedPolys.get(i).trigo;
		}
		return rots;
	}

	public double[][] getPackingPositions() {
		int size = packedPolys.size();
		double[][] ps = new double[size][];
		for (int i = 0; i < size; i++) {
			ps[i] = packedPolys.get(i).position;
		}
		return ps;
	}

	@Override
	public DalsooPack clone() {
		DalsooPack np = new DalsooPack(trigos, rotSteps, binWidth, binHeight, preferX);
		np.pendingPolys = pendingPolys;
		np.cntConvex = cntConvex;
		return np;
	}

	public boolean isEmpty() {
		return pendingPolys.isEmpty();
	}

}
