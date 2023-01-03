package whitegreen.dalsoo;

import java.util.ArrayList;

//Hao Hua, Southeast University, whitegreen@163.com

public class Strip implements Comparable<Strip> {
	
	public double[][] outps; // offset & add edge points , as referent point of placement, clock-wise
	public final double[][] inps; // original polygon, for 1. convex, 2. intersection
	public final double inarea;
	public double[] trigo; // cos & sin
	public double[] position;
	public final int id;

	public Strip(int id, double[][] original_poly, double spacing, Double segment_max_length) {
		this.id = id;
		if (spacing < 0) {
			throw new IllegalArgumentException("spacing must be >=0");
		}
		double area = M.area(original_poly);
		if (0 < area) {
			inps = M.clone(original_poly);
		} else {
			inps = new double[original_poly.length][];
			for (int i = 0; i < inps.length; i++) {
				inps[i] = original_poly[inps.length - 1 - i].clone();
			}
		}
		inarea = Math.abs(area);
		outps = M.offset(spacing, inps); // clockwise

		if (segment_max_length != null) {
			ArrayList<double[]> list = new ArrayList<>();
			for (int i = 0; i < outps.length; i++) {
				double[] pa = outps[i];
				double[] pb = outps[(i + 1) % outps.length];
				list.add(pa);
				double dis = M.dist(pa, pb);
				if (segment_max_length < dis) {
					int num = 1 + (int) (dis / segment_max_length);
					for (int j = 1; j < num; j++) {
						double s = (double) j / num;
						list.add(M.between(s, pa, pb));
					}
				}
			} // for
			outps = new double[list.size()][];
			for (int i = 0; i < outps.length; i++) {
				outps[i] = list.get(i);
			}
		} // if
	}

	@Override
	public int compareTo(Strip pl) {
		return Double.compare(this.inarea, pl.inarea);
	}

	public void fix_rotate_move(double[] cossin, double[] dv) { // finalize
		trigo = cossin.clone();
		position = dv.clone();
		double cos = trigo[0];
		double sin = trigo[1];
		for (int i = 0; i < inps.length; i++) {
			double[] p = inps[i];
			double x = cos * p[0] - sin * p[1];
			double y = sin * p[0] + cos * p[1];
			inps[i] = new double[] { dv[0] + x, dv[1] + y };
		}
		for (int i = 0; i < outps.length; i++) {
			double[] p = outps[i];
			double x = cos * p[0] - sin * p[1];
			double y = sin * p[0] + cos * p[1];
			outps[i] = new double[] { dv[0] + x, dv[1] + y };
		}
	}

}
