using NUnit.Framework;

namespace NetTopologySuite.IO.Oracle.Test
{
    /// <summary>
    /// 
    /// </summary>
    [TestFixture]
    public class OracleTest
    {

        private static readonly OracleGeometryReader or = new OracleGeometryReader(NtsGeometryServices.Instance);
        private static readonly WKTReader wr = new WKTReader { IsOldNtsCoordinateSyntaxAllowed = false };

        /// <summary>
        /// 
        /// </summary>
        [Test]
        public void CCWTestsOnPolygon()
        {
            string wrongCCW = "POLYGON((10 10, 10 20, 20 20, 20 10, 10 10),(5 5,6 5,6 6,5 6,5 5))";
            string correctCCW = "POLYGON((10 10, 20 10, 20 20, 10 20, 10 10),(5 5,5 6,6 6,6 5,5 5))";

            var geom1 = wr.Read(wrongCCW);
            var geom2 = wr.Read(correctCCW);

            var t = new OracleGeometryWriter().Write(geom1);
            var geom3 = or.Read(t);
            Assert.That(geom2.EqualsExact(geom3));
        }

        /// <summary>
        /// Tests all geometry types by transforming from wkt to oracle and back
        /// </summary>
        /// <param name="wkt"></param>
        /// <param name="srid"></param>
        [TestCase("POINT(10 10)", -1)]
        [TestCase("POINT(10 10)", 4326)]
        [TestCase("POINT Z(10 10 0)", -1)]
        [TestCase("POINT Z(10 10 20)", -1)]
        [TestCase("MULTIPOINT(11 12, 20 20)", -1)]
        [TestCase("MULTIPOINT Z(11 12 12, 20 20 20)", -1)]
        [TestCase("LINESTRING(10 10,20 20,50 50,34 34)", -1)]
        [TestCase("LINESTRING Z(10 10 20,20 20 20,50 50 50,34 34 34)", -1)]
        [TestCase("POLYGON((10 10,20 10,20 20,10 20,10 10))", -1)]
        [TestCase("POLYGON((10 10,20 10,20 20,10 20,10 10),(5 5,5 6,6 6,6 5,5 5))", -1)]
        [TestCase("POLYGON Z((10 10 0,20 10 0,20 20 0,10 20 0,10 10 0),(5 5 0,5 6 0,6 6 0,6 5 0,5 5 0))", -1)]
        [TestCase("MULTIPOLYGON(((10 10,20 10,20 20,20 10,10 10)),((10 10,20 10,20 20,20 10,10 10)))", -1)]
        [TestCase("MULTIPOLYGON(((10 10,20 10,20 20,10 20,10 10),(5 5,5 6,6 6,6 5,5 5)),((10 10,20 10,20 20,20 10,10 10)))", -1)]
        [TestCase("MULTIPOLYGON(((10 10,20 10,20 20,10 20,10 10),(5 5,5 6,6 6,6 5,5 5)),((10 10,20 10,20 20,20 10,10 10),(5 5,5 6,6 6,6 5,5 5)))", -1)]
        [TestCase("MULTIPOLYGON Z(((10 10 0,20 10 0,20 20 0,10 20 0,10 10 0),(5 5 0,5 6 0,6 6 0,6 5 0,5 5 0)),((10 10 0,20 10 0,20 20 0,20 10 0,10 10 0),(5 5 0,5 6 0,6 6 0,6 5 0,5 5 0)))", -1)]
        [TestCase("MULTILINESTRING((10 10,20 10,20 20,20 10),(5 5,5 6,6 6,6 5))", -1)]
        [TestCase("MULTILINESTRING((1 1, 2 1, 3 1), (1 2, 2 2, 3 2, 4 2), (1 3, 1 3, 3 3, 4 3))", -1)]
        [TestCase("MULTILINESTRING((1 1, 2 1, 3 1), (1 2, 2 2, 3 2, 4 2), (1 3, 1 3, 3 3, 4 3),(1 5, 2 5, 3 5),(1 6, 2 6, 3 6, 4 6))", -1)]
        [TestCase("MULTILINESTRING Z((10 10 5,20 10 5,20 20 0,20 10 0,10 10 0),(5 5 0,5 6 0,6 6 0,6 5 0,5 5 0))", -1)]
        [TestCase("GEOMETRYCOLLECTION(POLYGON((10 10,20 10,20 20,10 20,10 10)),POLYGON((30 10,40 10,40 20,30 20,30 10)))", -1)]
        [TestCase("GEOMETRYCOLLECTION(POLYGON((10 10,20 10,20 20,10 20,10 10),(5 5,5 6,6 6,6 5,5 5)))", -1)]
        [TestCase("GEOMETRYCOLLECTION(POLYGON((10 10,20 10,20 20,10 20,10 10),(5 5,5 6,6 6,6 5,5 5)),LINESTRING(10 10,20 20,50 50,34 34))", -1)]
        [TestCase("GEOMETRYCOLLECTION(POINT(10 10),LINESTRING(10 10,20 20,50 50,34 34))", -1)]
        [TestCase("GEOMETRYCOLLECTION(POINT(10 10),MULTIPOINT(11 12, 20 20))", -1)]
        public void BasicConversion(string wkt, int srid)
        {
            var geom = wr.Read(wkt);
            string parsed = geom.AsText();
            var regeom = wr.Read(parsed);
            string reparsed = regeom.AsText();

            geom.SRID = srid;
            regeom.SRID = srid;

            Assert.That(geom.EqualsExact(regeom));
            Assert.That(reparsed, Is.EqualTo(parsed));

            var t = new OracleGeometryWriter().Write(regeom);
            var regeom3 = or.Read(t);
            Assert.That(geom.EqualsExact(regeom3));
        }

        /// <summary>
        /// Tests geometry collection with multitypes 
        /// </summary>
        /// <param name="wkt"></param>
        /// <param name="wktresult"></param>
        /// <param name="srid"></param>
        [TestCase("GEOMETRYCOLLECTION(MULTIPOINT(11 12, 20 20))", "GEOMETRYCOLLECTION(MULTIPOINT(11 12, 20 20))", - 1)]
        [TestCase("GEOMETRYCOLLECTION(MULTIPOLYGON(((10 10,20 10,20 20,10 20,10 10),(5 5,5 6,6 6,6 5,5 5)),((10 10,20 10,20 20,20 10,10 10),(5 5,5 6,6 6,6 5,5 5))))", "GEOMETRYCOLLECTION(POLYGON((10 10,20 10,20 20,10 20,10 10),(5 5,5 6,6 6,6 5,5 5)),POLYGON((10 10,20 10,20 20,20 10,10 10),(5 5,5 6,6 6,6 5,5 5)))", - 1)]
        [TestCase("GEOMETRYCOLLECTION(MULTILINESTRING((10 10,20 10,20 20,10 20,10 10),(5 5,5 6,6 6,6 5,5 5)))", "GEOMETRYCOLLECTION(LINESTRING(10 10,20 10,20 20,10 20,10 10),LINESTRING(5 5,5 6,6 6,6 5,5 5))", -1)]
        public void CollectionConversion(string wkt, string wktresult, int srid)
        {
            var geom = wr.Read(wkt);
            var result = wr.Read(wktresult);

            geom.SRID = srid;
            result.SRID = srid;

            var t = new OracleGeometryWriter().Write(geom);
            var regeom = or.Read(t);
            Assert.That(result.EqualsExact(regeom));
        }

    }
}
