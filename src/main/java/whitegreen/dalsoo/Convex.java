package whitegreen.dalsoo;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Hao Hua, Southeast University, whitegreen@163.com
 */
class Convex {

	private static final double ZERO_AREA = 1E-6; // for mm
	
	List<double[]> convex; // convex is clockwise

	Convex() {
	}

	public Convex(PackedPoly strip) {
		double[][] points = strip.inpts;
		double area = strip.inarea;
		if (Math.abs(area) < ZERO_AREA) {
			throw new RuntimeException();
		}
		convex = new ArrayList<>();
		if (area > 0) {
			convex.add(points[0]);
			convex.add(points[1]);
			convex.add(points[2]);
		} else {
			convex.add(points[1]);
			convex.add(points[0]);
			convex.add(points[2]);
		}
		for (int i = 3; i < points.length; i++) {
			increment_hull(points[i]);
		}
	}

	@Override
	public Convex clone() {
		List<double[]> list = new ArrayList<>();
		list.addAll(convex);
		Convex con = new Convex();
		con.convex = list;
		return con;
	}

	public void increment_hull(double[] np) { // convex is clockwise
		if (MathUtil.inside(np, convex)) {
			return;
		}
		int size = convex.size();
		boolean[] visible = new boolean[size];
		boolean hasTure = false;
		for (int i = 0; i < size; i++) {
			double[] q0 = convex.get(i);
			double[] q1 = convex.get((i + 1) % size);
			double[] d = MathUtil.sub(q1, q0);
			double[] a = MathUtil.sub(np, q0);
			visible[i] = -d[1] * a[0] + d[0] * a[1] > 0; // n={-d[1], d[0]} on the left
			if (visible[i]) {
				hasTure = true;
			}
		}
		if (!hasTure) {
			return;
		}

		int T1st = -1;
		for (int i = 0; i < size; i++) { // find 1st visible
			int j = (i + 1) % size;
			if (!visible[i] && visible[j]) {
				T1st = j;
				break;
			}
		}
		if (0 > T1st) {
//			System.out.println(np[0] + "*" + np[1]);
//			for (double[] p : convex) {
//				System.out.println("new double[]{" + p[0] + "," + p[1] + "}");
//			}
//			for (boolean b : visible) {
//				System.out.println(b);
//			}
			throw new RuntimeException();
		}
		List<double[]> list = new ArrayList<>(size);
		list.add(convex.get(T1st));
		list.add(np);
		for (int i = 0; i < size; i++) {
			int j = (T1st + i + 1) % size;
			if (!visible[j]) {
				list.add(convex.get(j));
			}
		}
		convex = list;
	}

}
