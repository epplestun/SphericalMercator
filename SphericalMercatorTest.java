import com.buntplanet.buntbrain.persistence.leakfinder.util.SphericalMercator;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class SphericalMercatorTest {

  @Test
  public void spherical_mercator_bbbox_4326() {
    SphericalMercator sphericalMercator = new SphericalMercator();
    double[] bbox = sphericalMercator.bbox(1, 1, 1, false, "4326");

    assertEquals(0, bbox[0], 0.0);
    assertEquals(-85.05112877980659, bbox[1], 0.0);
    assertEquals(180, bbox[2], 0.0);
    assertEquals(0, bbox[3], 0.0);
  }

  @Test
  public void spherical_mercator_bbbox_900913() {
    SphericalMercator sphericalMercator = new SphericalMercator();
    double[] bbox = sphericalMercator.bbox(1, 1, 1, false, "900913");

    assertEquals(0, bbox[0], 0.0);
    assertEquals(-20037508.342789236, bbox[1], 0.0);
    assertEquals(20037508.342789244, bbox[2], 0.0);
    assertEquals(-7.081154551613622e-10, bbox[3], 0.0);
  }
}
