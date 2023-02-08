package whitegreen.dalsoo;

import java.util.ArrayList;
import java.util.List;

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
public class DalsooPack {

	List<Bin> bins;
	private final Bin initial;
	private final double binWidth, binHeight;

	/**
	 * Creates an packing. Call {@link #pack(boolean, boolean) packOneBin()} to
	 * start packing.
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
		initial = new Bin(polys, spacing, segment_max_length, rotSteps, width, height, hSkew);
		this.binWidth = width;
		this.binHeight = height;

		bins = new ArrayList<>();
	}

	/**
	 * Packs one bin with any remaining unpacked polygons. This operation does
	 * nothing if there are no more polygons to pack.
	 * 
	 * @param abey         whether to use the abey packing method
	 * @param largestFirst whether to pack polygons according to their area,
	 *                     largest-first
	 */
	public void packNextBin(boolean abey, boolean largestFirst) {
		Bin pack;
		if (bins.isEmpty()) {
			pack = initial;
		} else {
			pack = bins.get(bins.size() - 1).clone();
		}
		pack.pack(abey, largestFirst);
		bins.add(pack);
	}

	/**
	 * Packs all polygons, possibly over multiple bins.
	 * 
	 * @param abey         whether to use the abey packing method
	 * @param largestFirst whether to pack polygons according to their area,
	 *                     largest-first
	 */
	public void packAll(boolean abey, boolean largestFirst) {
		while (bins.isEmpty() || !bins.get(bins.size() - 1).isEmpty()) {
			packNextBin(abey, largestFirst);
		}
	}
	
	public List<double[][]> getPackedPolys(int maxHorizontalBins) {
		return getPackedPolys(maxHorizontalBins, 0);
	}

	public List<double[][]> getPackedPolys(int maxHorizontalBins, double binSpacing) {
		List<double[][]> polys = new ArrayList<>();
		for (int i = 0; i < bins.size(); i++) {
			int column = i % maxHorizontalBins;
			int row = i / maxHorizontalBins;
			double xOffset = column * (binWidth + binSpacing);
			double yOffset = row * (binHeight + binSpacing);
			final double[] vector = { xOffset, yOffset };
			bins.get(i).getPackedPolygons().forEach(p -> {
				polys.add(MathUtil.translate(vector, p.getPackedCoordinates()));
			});
		}
		return polys;
	}

	public List<Bin> getPackedBins() {
		return bins;
	}

}
