import org.apache.commons.lang.ArrayUtils;

import java.util.HashMap;
import java.util.Map;

public class SphericalMercator {
  private final double EPSLN = 1.0e-10;
  private final double D2R = Math.PI / 180;
  private final double R2D = 180 / Math.PI;
  private final double A = 6378137.0;
  private final double MAXEXTENT = 20037508.342789244;
  private final Map<Double, HashMap<String, double[]>> cache = new HashMap<>();
  private final double[] Bc;
  private final double[] Cc;
  private final double[] zc;
  private final double[] Ac;

  private double size;

  public SphericalMercator() {
    this.size = 256;
    if (this.cache.get(this.size) == null) {
      double size = this.size;

      HashMap<String, double[]> c = new HashMap<String, double[]>() {{
        put("Bc", new double[30]);
        put("Cc", new double[30]);
        put("zc", new double[30]);
        put("Ac", new double[30]);
      }};

      this.cache.put(size, c);

      double[] Bc = c.get("Bc");
      double[] Cc = c.get("Cc");
      double[] zc = c.get("zc");
      double[] Ac = c.get("Ac");

      for (int d = 0; d < 30; d++) {
        Bc[d] = size / 360;
        Cc[d] = size / (2 * Math.PI);
        zc[d] = size / 2;
        Ac[d] = size;

        size *= 2;
      }
    }

    this.Bc = this.cache.get(this.size).get("Bc");
    this.Cc = this.cache.get(this.size).get("Cc");
    this.zc = this.cache.get(this.size).get("zc");
    this.Ac = this.cache.get(this.size).get("Ac");
  }

  public double[] px(double[] ll, int zoom) {
    double d = this.zc[zoom];
    double f = Math.min(Math.max(Math.sin(D2R * ll[1]), -0.9999), 0.9999);
    double x = Math.round(d + ll[0] * this.Bc[zoom]);
    double y = Math.round(d + 0.5 * Math.log((1 + f) / (1 - f)) * (-this.Cc[zoom]));

    double[] px = new double[2];
    px[0] = x;
    px[1] = y;

    return px;
  }

  public double[] ll(double[] px, int zoom) {
    double g = (px[1] - this.zc[zoom]) / (-this.Cc[zoom]);
    double lon = (px[0] - this.zc[zoom]) / this.Bc[zoom];
    double lat = R2D * (2 * Math.atan(Math.exp(g)) - 0.5 * Math.PI);

    double[] ll = new double[2];
    ll[0] = lon;
    ll[1] = lat;

    return ll;
  }

  public double[] bbox(double x, double y, int zoom, boolean tms_style, String srs) {
    // Convert xyz into bbox with srs WGS84
    if (tms_style) {
      y = (Math.pow(2, zoom) - 1) - y;
    }

    double[] ll = new double[2];
    ll[0] = x * this.size;
    ll[1] = (+y + 1) * this.size;

    double[] ur = new double[2];
    ur[0] = (+x + 1) * this.size;
    ur[1] = y * this.size;

    double[] bbox = ArrayUtils.addAll(this.ll(ll, zoom), this.ll(ur, zoom));

    // If web mercator requested reproject to 900913.
    if (srs.equals("900913")) {
      return this.convert(bbox, "900913");
    } else {
      return bbox;
    }
  }

  public double[] xyz(double[] bbox, int zoom, boolean tms_style, String srs) {
    // If web mercator provided reproject to WGS84.
    if (srs.equals("900913")) {
      bbox = this.convert(bbox, "WGS84");
    }

    double[] ll = new double[2];
    ll[0] = bbox[0];
    ll[1] = bbox[1];

    double[] ur = new double[2];
    ur[0] = bbox[2];
    ur[1] = bbox[3];

    double[] px_ll = this.px(ll, zoom);
    double[] px_ur = this.px(ur, zoom);

    double[] x = new double[2];
    x[0] = Math.floor(px_ll[0] / this.size);
    x[1] = Math.floor((px_ur[0] - 1) / this.size);

    double[] y = new double[2];
    y[0] = Math.floor(px_ur[1] / this.size);
    y[1] = Math.floor((px_ll[1] - 1) / this.size);

    double[] bounds = new double[4];
    bounds[0] = Math.min(x[0], x[1]) < 0 ? 0 : Math.min(x[0], x[1]); // minX
    bounds[1] = Math.min(y[0], y[1]) < 0 ? 0 : Math.min(y[0], y[1]); // minY
    bounds[2] = Math.max(x[0], x[1]);                                // maxX
    bounds[3] = Math.max(y[0], y[1]);                                // maxY

    if (tms_style) {
      bounds[1] = (Math.pow(2, zoom) - 1) - bounds[3];
      bounds[3] = (Math.pow(2, zoom) - 1) - bounds[1];
    }

    return bounds;
  }

  public double[] convert(double[] bbox, String to) {
    double[] a = bbox.clone();
    double[] c = new double[2];
    c[0] = a[0];
    c[1] = a[1];
    double[] d = new double[2];
    d[0] = a[2];
    d[1] = a[3];

    if (to.equals("900913")) {
      return ArrayUtils.addAll(this.forward(c), this.forward(d));
    } else {
      return ArrayUtils.addAll(this.inverse(c), this.inverse(d));
    }
  }

  public double[] forward(double[] ll) {
    double[] xy = new double[2];
    xy[0] = A * ll[0] * D2R;
    xy[1] = A * Math.log(Math.tan((Math.PI * 0.25) + (0.5 * ll[1] * D2R)));

    return xy;
  }

  public double[] inverse(double[] xy) {
    double[] result = new double[2];
    result[0] = (xy[0] * R2D / A);
    result[1] = ((Math.PI * 0.5) - 2.0 * Math.atan(Math.exp(-xy[1] / A))) * R2D;

    return result;
  }
}
