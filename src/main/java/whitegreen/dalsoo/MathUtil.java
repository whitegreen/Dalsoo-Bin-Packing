package whitegreen.dalsoo;

import java.util.ArrayList;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.buffer.BufferOp;
import org.locationtech.jts.operation.buffer.BufferParameters;

import net.jafama.FastMath;

//Hao Hua, Southeast University, whitegreen@163.com

public final class MathUtil {

	private static final double denominator_lim = 1e-10;
	private static final BufferParameters bp = new BufferParameters(8, BufferParameters.CAP_FLAT, BufferParameters.JOIN_MITRE,
			BufferParameters.DEFAULT_MITRE_LIMIT);

	// math
	// *************************************************************************************
	public static int powerInt(int a) {
		int v = 1;
		for (int i = 2; i <= a; i++) {
			v *= i;
		}
		return v;
	}

	// vector **********************************************************************
	public static void assign(double[] a, double[] b) {
		for (int i = 0; i < a.length; i++) {
			a[i] = b[i];
		}
	}

	public static double[] mean(double[][] ps) {
		double[] re = new double[ps[0].length];
		for (double[] element : ps) {
			for (int j = 0; j < re.length; j++) {
				re[j] += element[j];
			}
		}
		for (int j = 0; j < re.length; j++) {
			re[j] /= ps.length;
		}
		return re;
	}

	public static double dist_sq(double[] a, double[] b, int len) {
		double sum = 0;
		for (int i = 0; i < len; i++) {
			double d = a[i] - b[i];
			sum += d * d;
		}
		return sum;
	}

	public static double dist_sq(double[] a, double[] b) {
		double sum = 0;
		for (int i = 0; i < a.length; i++) {
			double d = a[i] - b[i];
			sum += d * d;
		}
		return sum;
	}

	public static double dist(double[] a, double[] b) {
		double sum = 0;
		for (int i = 0; i < a.length; i++) {
			double d = a[i] - b[i];
			sum += d * d;
		}
		return FastMath.sqrt(sum);
	}

	public static double triArea(double[] a, double[] b) {
		double[] cross = MathUtil.cross(a, b);
		return MathUtil.mag(cross) * 0.5;
	}

	public static double mag(double[] a) {
		double sum = 0;
		for (double element : a) {
			sum += element * element;
		}
		return FastMath.sqrt(sum);
	}

	public static double mag_sq(double[] a) {
		double sum = 0;
		for (double element : a) {
			sum += element * element;
		}
		return sum;
	}

	public static void swap(double[] a, double[] b) {
		double[] na = a.clone();
		double[] nb = b.clone();
		for (int i = 0; i < a.length; i++) {
			a[i] = nb[i];
			b[i] = na[i];
		}
	}

	public static double[][] clone(double[][] vs) {
		double[][] ps = new double[vs.length][];
		for (int i = 0; i < vs.length; i++) {
			ps[i] = vs[i].clone();
		}
		return ps;
	}

	public static double[] sub(double[] a, double[] b) {
		double[] p = new double[a.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = a[i] - b[i];
		}
		return p;
	}

	public static double[] sub(double[] a, double[] b, int len) {
		double[] p = new double[a.length];
		for (int i = 0; i < len; i++) {
			p[i] = a[i] - b[i];
		}
		return p;
	}

	public static void _sub(double[] a, double[] b) {
		for (int i = 0; i < a.length; i++) {
			a[i] -= b[i];
		}
	}

	public static double[] add(double[] a, double[] b) {
		double[] p = new double[a.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = a[i] + b[i];
		}
		return p;
	}

	public static double[] add(double[][] vs) {
		int d = vs[0].length;
		double[] p = new double[d];
		for (double[] element : vs) {
			for (int j = 0; j < d; j++) {
				p[j] += element[j];
			}
		}
		return p;
	}

	public static double[] add(ArrayList<double[]> vs) {
		int d = vs.get(0).length;
		double[] p = new double[d];
		for (double[] v : vs) {
			for (int j = 0; j < d; j++) {
				p[j] += v[j];
			}
		}
		return p;
	}

	public static double[] add(double sa, double[] a, double sb, double[] b) {
		double[] p = new double[a.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = sa * a[i] + sb * b[i];
		}
		return p;
	}

	public static void _add(double[] a, double[] b) {
		for (int i = 0; i < a.length; i++) {
			a[i] += b[i];
		}
	}

	public static double[] scale(double s, double[] a) {
		double[] p = new double[a.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = a[i] * s;
		}
		return p;
	}

	public static void _scale(double s, double[] a) {
		for (int i = 0; i < a.length; i++) {
			a[i] *= s;
		}
	}

	public static double[] scaleTo(double scale, double[] a) {
		double sc = scale / mag(a);
		double[] p = new double[a.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = a[i] * sc;
		}
		return p;
	}

	public static void _scaleTo(double scale, double[] a) {
		double sc = scale / mag(a);
		for (int i = 0; i < a.length; i++) {
			a[i] *= sc;
		}
	}

	public static double[] normalize(double[] v) {
		double mag = mag(v);
		double[] p = new double[v.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = v[i] / mag;
		}
		return p;
	}

	public static void _normalize(double[] v) {
		double mag = mag(v);
		for (int i = 0; i < v.length; i++) {
			v[i] /= mag;
		}
	}

	public static double[] between(double s, double[] a, double[] b) {
		double[] p = new double[a.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = (1 - s) * a[i] + s * b[i];
		}
		return p;
	}

	public static double dot(double[] a, double[] b) {
		double sum = 0;
		for (int i = 0; i < a.length; i++) {
			sum += a[i] * b[i];
		}
		return sum;
	}

	public static double[] cross(double[] a, double[] b) {
		double x = a[1] * b[2] - a[2] * b[1];
		double y = a[2] * b[0] - a[0] * b[2];
		double z = a[0] * b[1] - a[1] * b[0];
		double[] v = { x, y, z };
		return v;
	}

	static double[] center(double[][] ps) {
		int len = ps[0].length;
		double[] sum = new double[len];
		for (double[] element : ps) {
			for (int j = 0; j < len; j++) {
				sum[j] += element[j];
			}
		}
		for (int j = 0; j < len; j++) {
			sum[j] /= ps.length;
		}
		return sum;
	}

	// linear algebra matrix
	// **********************************************************************
	public static double[][] transpose(double[][] v) {
		int row = v.length;
		int col = v[0].length;
		double[][] re = new double[col][row];
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				re[j][i] = v[i][j];
			}
		}
		return re;
	}

	public static double[][] sub(double[][] a, double[][] b) {
		int row = a.length;
		int col = a[0].length;
		double[][] re = new double[row][col];
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				re[i][j] = a[i][j] - b[i][j];
			}
		}
		return re;
	}

	public static double[][] add(double[][] a, double[][] b) {
		int row = a.length;
		int col = a[0].length;
		double[][] re = new double[row][col];
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				re[i][j] = a[i][j] + b[i][j];
			}
		}
		return re;
	}

//	public static double[][] mul(double[][] a, double[][] b) {
//		int row = a.length;
//		int col = a[0].length;
//		double[][] re = new double[row][col];
//		for (int i = 0; i < row; i++) {
//			for (int j = 0; j < col; j++) {
//				for (int k = 0; k < col; k++)
//					re[i][j] += a[i][k] * b[k][j];
//			}
//		}
//		return re;
//	}
	public static double[][] mul(double[][] a, double[][] b) {
		int row = a.length;
		int col = b[0].length;
		double[][] re = new double[row][col];
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				for (int k = 0; k < b.length; k++) {
					re[i][j] += a[i][k] * b[k][j];
				}
			}
		}
		return re;
	}

	public static double[][] scale(double[][] a, double s) {
		int row = a.length;
		int col = a[0].length;
		double[][] re = new double[row][col];
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				re[i][j] = s * a[i][j];
			}
		}
		return re;
	}

	public static double[][] mul(double[][] m, int iteration) {
		double[][] a = m;
		for (int i = 0; i < iteration - 1; i++) {
			a = MathUtil.mul(a, m);
		}
		return a;
	}

	public static double trace(double[][] m) {
		double sum = 0;
		for (int i = 0; i < m.length; i++) {
			sum += m[i][i];
		}
		return sum;
	}

	public static double[] mul(double[][] m, double[] p) {
		double[] v = new double[p.length];
		for (int i = 0; i < p.length; i++) { // row
			for (int j = 0; j < p.length; j++) { // col
				v[i] += m[i][j] * p[j];
			}
		}
		return v;
	}

	public static double[][] I(int d) {
		double[][] m = new double[d][d];
		for (int i = 0; i < d; i++) {
			m[i][i] = 1;
		}
		return m;
	}

	public static double[][] rotate(double[] trigo, double[][] ps) {
		double[][] arr = new double[ps.length][];
		for (int i = 0; i < ps.length; i++) {
			arr[i] = rotate(trigo, ps[i]);
		}
		return arr;
	}

	public static double[] rotate(double[] trigo, double[] p) {
		double cos = trigo[0];
		double sin = trigo[1];
		double x = cos * p[0] - sin * p[1];
		double y = sin * p[0] + cos * p[1];
		return new double[] { x, y };
	}

	public static double[][] move(double[] v, double[][] ps) {
		double[][] arr = new double[ps.length][];
		for (int i = 0; i < ps.length; i++) {
			arr[i] = add(ps[i], v);
		}
		return arr;
	}

	// computational geometry
	// ************************************************************************
	public static boolean intersect_boundBox(double[][] poly1, double[][] poly2) {
		double[] a = boundBox(poly1);
		double[] b = boundBox(poly2);
		if (a[1] < b[0] || a[0] > b[1] || a[3] < b[2] || a[2] > b[3]) {
			return false;
		}
		return true;
	}

	/**
	 * 
	 * @param a [minx, maxx, miny, maxy]
	 * @param b [minx, maxx, miny, maxy]
	 * @return
	 */
	public static boolean intersect_boundBoxes(double[] a, double[] b) {
		if (a[1] < b[0] || a[0] > b[1] || a[3] < b[2] || a[2] > b[3]) {
			return false;
		}
		return true;
	}

	public static double[] boundBox(double[][] poly) { // orthogonal, minx, max, miny, maxy
		double[] v = poly[0];
		double minx = v[0];
		double maxx = v[0];
		double miny = v[1];
		double maxy = v[1];
		for (double[] p : poly) {
			if (minx > p[0]) {
				minx = p[0];
			}
			if (maxx < p[0]) {
				maxx = p[0];
			}

			if (miny > p[1]) {
				miny = p[1];
			}
			if (maxy < p[1]) {
				maxy = p[1];
			}
		}
		return new double[] { minx, maxx, miny, maxy };
	}

	public static double[][] buffer(final double[][] poly, final double distance) {
		double[][] ps = new double[poly.length][];
		for (int i = 0; i < poly.length; i++) {
			double[] p = poly[i];
			double[] va = sub(poly[(i + 1) % poly.length], p);
			double[] vb = sub(p, poly[(i - 1 + poly.length) % poly.length]);
			double[] na = { -va[1], va[0] };
			double[] nb = { -vb[1], vb[0] };
			double[] pa = distance == 0 ? p : add(p, scaleTo(distance, na));
			double[] pb = distance == 0 ? p : add(p, scaleTo(distance, nb));
			ps[i] = lineIntersect(pa, va, pb, vb);
			if (ps[i] == null) {
				/*
				 * Fall back to JTS buffer if lineIntersect returns null for any intersection
				 * (occurs on rather concave shapes).
				 */
				return toArray(BufferOp.bufferOp(toPolygon(poly), distance, bp).getCoordinates());
			}
		}
		return ps;
	}
	
	private static final GeometryFactory GEOM_FACT = new GeometryFactory();
	
	private static Polygon toPolygon(double[][] poly) {
		Coordinate[] coords = new Coordinate[poly.length + 1];
		for (int i = 0; i < coords.length - 1; i++) {
			coords[i] = new Coordinate(poly[i][0], poly[i][1]);
		}
		coords[coords.length - 1] = new Coordinate(poly[0][0], poly[0][1]); // close
		return GEOM_FACT.createPolygon(coords);
	}

	private static double[][] toArray(Coordinate[] coords) {
		double[][] poly = new double[coords.length - 1][2];
		for (int i = 0; i < coords.length - 1; i++) { // -1, unclose
			poly[i][0] = coords[i].x;
			poly[i][1] = coords[i].y;
		}
		return poly;
	}

	public static double[] point_Project_Plane(double[] q, double[] p, double[] n) { // q is the point
		double d = -dot(p, n);
		double[] tmp = scale(dot(q, n) + d, n);
		return sub(q, tmp);
	}

	public static double[] vector_Project_Plane(double[] v, double[] p, double[] n) { // v is the vector
		double[] tmp = scale(dot(v, n), n);
		return sub(v, tmp);
	}

	public static double[][] circleLine_intersection(double[] cnt, double r, double[] p, double[] d) { // two intersections or
																										// nothing
		// d must be normalized
		System.out.println(mag_sq(d));
		double[] m = sub(p, cnt);
		double b = 2 * dot(m, d);
		double c = mag_sq(m) - r * r;
		double b4c = b * b - 4 * c;
		if (b4c < 0) {
			return null;
		}
		double delta = FastMath.sqrt(b4c);
		double s1 = 0.5 * (-b + delta);
		double s2 = 0.5 * (-b - delta);
		double[][] re = { add(p, scale(s1, d)), add(p, scale(s2, d)) };
		return re;
	}

	public static double[] planes_Intersect(double[] r1, double[] n1, double[] r2, double[] n2) { // wiki plane
		double h1 = dot(n1, r1);
		double h2 = dot(n2, r2);
		double dn12 = dot(n1, n2);
		double fem = 1 - dn12 * dn12;
		if (FastMath.abs(fem) < denominator_lim) {
			// return M.between(0.5, r1, r2);
			throw new RuntimeException();
		}
		double c1 = (h1 - h2 * dn12) / fem;
		double c2 = (h2 - h1 * dn12) / fem;
		return add(c1, n1, c2, n2);
		// double[] direction = crossProduct(n1, n2);
		// _normalize(direction);
	}

	public static double[] planeLine_Intersect(double[] n, double[] point, double[] p, double[] dir) { // n={a b c}
		return planeLine_Intersect(n, -dot(n, point), p, dir);
	}

	public static double[] planeLine_Intersect(double[] n, double d, double[] p, double[] dir) { // n={a b c} n&d: plane, p&dir:
																									// line
		double fz = (dot(n, p)) + d;
		double fm = dot(n, dir);
		if (FastMath.abs(fm) < denominator_lim) {
			return null;
		}
		double t = -fz / fm;
		return add(p, scale(t, dir));
	}

	public static double[] lineIntersect3D(double[] p0, double[] d0, double[] p1, double[] d1) { // project p1d1 to p0d0, d0, d1
																									// normalized
		if (FastMath.abs(dot(d0, d1)) < denominator_lim) { // parallel
			return null;
		}
		double[] dx = d0; // mag=1
		double[] dz = cross(d0, d1);
		double[] dy = cross(dx, dz);
		_normalize(dy); // mag=1;

		double[] arrow = sub(p1, p0);
		double x = dot(dx, arrow);
		double y = dot(dy, arrow);
		double vx = dot(dx, d1);
		double vy = dot(dy, d1);
		double xx = x - (y * vx) / vy;
		return add(p0, scale(xx, dx));
	}

	public static double[] lineIntersect(double[] p0, double[] n0, double[] p1, double[] n1) { // in 2d
		double cross_base = kross(n0, n1);
		if (FastMath.abs(cross_base) < denominator_lim) { // parallel
			return null;
		}
		double[] d = sub(p1, p0);
		double s = kross(d, n1) / cross_base;
		return new double[] { p0[0] + s * n0[0], p0[1] + s * n0[1] };
	}

	private static double kross(double[] a, double[] b) {
		return a[0] * b[1] - a[1] * b[0];
	}

	public static double area(final double[][] ps) { // signed
		double sum = 0;
		for (int i = 0; i < ps.length; i++) {
			double[] pa = ps[i];
			double[] pb = ps[(i + 1) % ps.length];
			sum += pa[1] * pb[0] - pa[0] * pb[1];
		}
		return 0.5f * sum;
	}

	public static double areaAbs(ArrayList<double[]> ps) { // signed
		double sum = 0;
		int size = ps.size();
		for (int i = 0; i < size; i++) {
			double[] pa = ps.get(i);
			double[] pb = ps.get((i + 1) % size);
			sum += pa[1] * pb[0] - pa[0] * pb[1];
		}
		return FastMath.abs(0.5f * sum);
	}

	public static double areaAbs(double[][] ps) {
		return FastMath.abs(area(ps));
	}

	public static double[][] getRotateMatrix(double[] xyz) { // double rx, double ry, double rz
		double a = -xyz[2]; // rz
		double b = -xyz[1];
		double c = -xyz[0];
		double ca = FastMath.cos(a);
		double sa = FastMath.sin(a);
		double cb = FastMath.cos(b);
		double sb = FastMath.sin(b);
		double cc = FastMath.cos(c);
		double sc = FastMath.sin(c);
		double[][] tt = new double[3][];
		tt[0] = new double[] { ca * cb, sa * cc + ca * sb * sc, sa * sc - ca * sb * cc };
		tt[1] = new double[] { -sa * cb, ca * cc - sa * sb * sc, ca * sc + sa * sb * cc };
		tt[2] = new double[] { sb, -cb * sc, cb * cc };
		return tt;
	}

	public static double[] rotateX(double theta, double[] p) {
		double c = FastMath.cos(theta);
		double s = FastMath.sin(theta);
		double[][] m = { { 1, 0, 0 }, { 0, c, -s }, { 0, s, c } };
		return mul(m, p);
	}

	public static double[] rotateY(double theta, double[] p) {
		double c = FastMath.cos(theta);
		double s = FastMath.sin(theta);
		double[][] m = { { c, 0, s }, { 0, 1, 0 }, { -s, 0, c } };
		return mul(m, p);
	}

	public static double[] rotateZ(double theta, double[] p) {
		double c = FastMath.cos(theta);
		double s = FastMath.sin(theta);
		double[][] m = { { c, -s, 0 }, { s, c, 0 }, { 0, 0, 1 } };
		return mul(m, p);
	}

	public static double[] rotate_Rodriguez(double[] v, double theta, double[] dir) {
		double cos = FastMath.cos(theta);
		double sin = FastMath.sin(theta);
		double[] v1 = MathUtil.scale(cos, v);
		double[] v2 = MathUtil.scale((1 - cos) * dot(v, dir), dir);
		double[] v3 = MathUtil.scale(sin, cross(dir, v));
		return MathUtil.add(v1, MathUtil.add(v2, v3));
	}

	public static boolean inside(double[] p, double[][] vs) {
		int i, j = vs.length - 1;
		boolean oddNodes = false;
		for (i = 0; i < vs.length; i++) {
			if ((vs[i][1] < p[1] && vs[j][1] >= p[1] || vs[j][1] < p[1] && vs[i][1] >= p[1])
					&& (vs[i][0] <= p[0] || vs[j][0] <= p[0])) {
				if (vs[i][0] + (p[1] - vs[i][1]) / (vs[j][1] - vs[i][1]) * (vs[j][0] - vs[i][0]) < p[0]) {
					oddNodes = !oddNodes;
				}
			}
			j = i;
		}
		return oddNodes;
	}

	public static boolean inside(double[] p, ArrayList<double[]> vs) {
		int size = vs.size();
		int i, j = size - 1;
		boolean oddNodes = false;
		for (i = 0; i < size; i++) {
			double[] vi = vs.get(i);
			double[] vj = vs.get(j);
			if ((vi[1] < p[1] && vj[1] >= p[1] || vj[1] < p[1] && vi[1] >= p[1]) && (vi[0] <= p[0] || vj[0] <= p[0])) {
				if (vi[0] + (p[1] - vi[1]) / (vj[1] - vi[1]) * (vj[0] - vi[0]) < p[0]) {
					oddNodes = !oddNodes;
				}
			}
			j = i;
		}
		return oddNodes;
	}

	/**
	 * Determines whether two line segments intersect.
	 * 
	 * @param a line 1 p1
	 * @param b line 1 p2
	 * @param c line 2 p1
	 * @param d line 2 p2
	 * @return
	 */
	public static boolean intersection(final double[] a, final double[] b, final double[] c, final double[] d) {
		double[][] abd = { a, b, d };
		double[][] abc = { a, b, c };
		double[][] cda = { c, d, a };
		double[][] cdb = { c, d, b };
		return area(abd) * area(abc) < 0 && area(cda) * area(cdb) < 0;
	}

	public static boolean intersectionFast(final double[] a, final double[] b, final double[] c, final double[] d) {
		return triangleArea(a, b, d) * triangleArea(a, b, c) < 0 && triangleArea(c, d, a) * triangleArea(c, d, b) < 0;
	}

	public static double triangleArea(final double[] a, final double[] b, final double c[]) {
		return (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1])); // omit /2
	}

	// special
	// *************************************************************************************

	public static double[][] exp(double[][] m, int iteration) {
		double[][] re = I(m.length);
		for (int i = 0; i < iteration; i++) {
			double[][] tmp = MathUtil.scale(mul(m, i + 1), 1.0 / powerInt(i + 1));
			re = MathUtil.add(re, tmp);
		}
		return re;
	}

	public static double[][] ln(double[][] rotate_matrix) {// group to algebra
		double tr = 0;
		for (int i = 0; i < rotate_matrix.length; i++) {
			tr += rotate_matrix[i][i];
		}
		double theta = FastMath.acos(0.5 * (tr - 1));
		double[][] t = MathUtil.sub(rotate_matrix, MathUtil.transpose(rotate_matrix));
		double[][] re = MathUtil.scale(t, 0.5 * theta / FastMath.sin(theta)); // sin(0)*************
		return re;
	}

	public static double[] algebra_ele(double[][] mat) {
		double[] ele = { -mat[1][2], mat[0][2], -mat[0][1] };
		return ele;
	}

	public static double[][] algebra_matrix(double[] ele) {
		double[][] mat = { { 0, -ele[2], ele[1] }, { ele[2], 0, -ele[0] }, { -ele[1], ele[0], 0 } };
		return mat;
	}

	public static double[][] exp(double[] ele) { // algebra to group
		double sq_theta = MathUtil.mag_sq(ele);
		double theta = FastMath.sqrt(sq_theta);
		double[][] ele_mat = algebra_matrix(ele);
		double[][] re = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
		re = MathUtil.add(re, MathUtil.scale(ele_mat, FastMath.sin(theta) / theta));
		return MathUtil.add(re, MathUtil.scale(MathUtil.mul(ele_mat, ele_mat), (1 - FastMath.cos(theta)) / sq_theta));
	}

}
