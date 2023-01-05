package whitegreen.dalsoo;

import java.util.ArrayList;

/**
 * 
 * @author Hao Hua, Southeast University, whitegreen@163.com
 *
 */
public class PackedPoly implements Comparable<PackedPoly> {

	public double[][] outpts; // offset & add edge points , as referent point of placement, clock-wise
	/** Coordinates of this poly's original polygon */
	public final double[][] inpts; // original polygon, for 1. convex, 2. intersection
	public final double inarea;
	public double[] trigo; // cos & sin
	public double[] position;
	public final int id;

	PackedPoly(int id, double[][] original_poly, double spacing, Double segmentMaxLength) {
		this.id = id;
		if (spacing < 0) {
			throw new IllegalArgumentException("spacing must be >=0");
		}
		double area = MathUtil.area(original_poly);
		if (0 < area) {
			inpts = MathUtil.clone(original_poly);
		} else {
			inpts = new double[original_poly.length][];
			for (int i = 0; i < inpts.length; i++) {
				inpts[i] = original_poly[inpts.length - 1 - i].clone();
			}
		}

		inarea = Math.abs(area);
		outpts = MathUtil.buffer(inpts, spacing); // clockwise

		if (segmentMaxLength != null && segmentMaxLength > 0) {
			ArrayList<double[]> list = new ArrayList<>();
			for (int i = 0; i < outpts.length; i++) {
				double[] pa = outpts[i];
				double[] pb = outpts[(i + 1) % outpts.length];
				list.add(pa);
				double dis = MathUtil.dist(pa, pb);
				if (segmentMaxLength < dis) {
					int num = 1 + (int) (dis / segmentMaxLength);
					for (int j = 1; j < num; j++) {
						double s = (double) j / num;
						list.add(MathUtil.between(s, pa, pb));
					}
				}
			} // for
			outpts = new double[list.size()][];
			for (int i = 0; i < outpts.length; i++) {
				outpts[i] = list.get(i);
			}
		} // if
	}

	void fix_rotate_move(double[] cossin, double[] dv) { // finalize
		trigo = cossin.clone();
		position = dv.clone();
		double cos = trigo[0];
		double sin = trigo[1];
		for (int i = 0; i < inpts.length; i++) {
			double[] p = inpts[i];
			double x = cos * p[0] - sin * p[1];
			double y = sin * p[0] + cos * p[1];
			inpts[i] = new double[] { dv[0] + x, dv[1] + y };
		}
		for (int i = 0; i < outpts.length; i++) {
			double[] p = outpts[i];
			double x = cos * p[0] - sin * p[1];
			double y = sin * p[0] + cos * p[1];
			outpts[i] = new double[] { dv[0] + x, dv[1] + y };
		}
	}

	double[] bb;
	double[] centroid;

	/**
	 * Call when placed in the bin. Will pre-compute and save values used for
	 * overlap detection against future polygons.
	 */
	void place() {
		bb = MathUtil.boundBox(inpts);
		centroid = MathUtil.center(inpts);
	}

	@Override
	public int compareTo(PackedPoly pl) {
		return Double.compare(this.inarea, pl.inarea);
	}

}
