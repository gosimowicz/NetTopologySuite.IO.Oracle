using System;
using System.Collections.Generic;
using System.Diagnostics;
using NetTopologySuite.Geometries;
using NetTopologySuite.IO.Sdo;

namespace NetTopologySuite.IO
{
    /// <summary>
    /// A class for reading Oracle <see cref="SdoGeometry"/> objects
    /// </summary>
    public class OracleGeometryReader
    {
        private const int NullDimension = -1;

        /// <summary>
        /// Creates an <c>OracleGeometryReader</c> using <see cref="GeometryFactory.Default"/>
        /// </summary>
        [Obsolete("Use constructor with NtsGeometryServices argument")]
        public OracleGeometryReader()
            : this(GeometryFactory.Default)
        { }

        /// <summary>
        /// Creates an <c>OracleGeometryReader</c> using the provided <see cref="GeometryFactory"/> object.
        /// </summary>
        /// <param name="factory">The geometry factory to use.</param>
        /// <remarks>
        /// Note that this factory will only be used as a template to create a new
        /// <see cref="NtsGeometryServices"/> object which is then used internally.
        /// </remarks>
        [Obsolete("Use constructor with NtsGeometryServices argument")]
        public OracleGeometryReader(GeometryFactory factory)
            : this(new NtsGeometryServices(factory.CoordinateSequenceFactory, factory.PrecisionModel, factory.SRID))
        {
        }

        /// <summary>
        /// Creates an <c>OracleGeometryReader</c> using the provided <see cref="NtsGeometryServices"/> object.
        /// </summary>
        public OracleGeometryReader(NtsGeometryServices services)
        {
            _services = services ?? NtsGeometryServices.Instance;
        }

        private readonly NtsGeometryServices _services;

        /// <summary>
        /// Gets or sets a value indicating the number of dimensions that are used
        /// </summary>
        public int Dimension { get; set; } = -1;

        /// <summary>
        /// Tranlates an Oracle <see cref="SdoGeometry"/> to a <c>Geometry</c>.
        /// </summary>
        /// <param name="geom">The Oracle geometry</param>
        /// <returns>The translated geometry</returns>
        public Geometry Read(SdoGeometry geom)
        {

            //Note: Returning null for null Datum
            if (geom == null)
                return null;

            Debug.Assert(geom.SdoGtype.HasValue);
            int gType = (int)geom.SdoGtype;

            Debug.Assert(geom.Sdo_Srid.HasValue);
            int srid = (int)geom.Sdo_Srid;

            var point = geom.Point;
            var factory = _services.CreateGeometryFactory(srid);

            var retVal = Create(factory, gType, point, geom.ElemArray, geom.OrdinatesArray);

            if (retVal == null)
            {
                return null;
            }
            retVal.SRID = srid;

            return retVal;
        }

        private Geometry Create(GeometryFactory factory, int gType, SdoPoint point, double[] elemInfo, double[] ordinates)
        {
            int lrs = (gType % 1000) / 100;

            // find the dimension: represented by the smaller of the two dimensions
            int dim;
            if (Dimension != NullDimension)
            {
                dim = Dimension;
            }
            else
            {
                dim = Math.Min(gType / 1000, 3);
            }

            if (dim == 0)
                return null;

            if (dim < 2)
            {
                throw new ArgumentException("Dimension D:" + dim + " is not valid for JTS. " +
                                            "Either specify a dimension or use Oracle Locator Version 9i or later");
            }

            // extract the geometry template type
            // this is represented as the rightmost two digits
            int geomTemplate = gType - (dim * 1000) - (lrs * 100);

            //CoordinateSequence coords = null;
            List<Coordinate> coords;

            if (lrs == 0 && geomTemplate == 1 && point != null && elemInfo == null)
            {
                // Single Coordinate Type Optimization
                Debug.Assert(point.X != null, "point.X != null");
                Debug.Assert(point.Y != null, "point.Y != null");
                if (dim == 2)
                {
                    coords = Coordinates(dim, lrs, geomTemplate, new[] { point.X.Value, point.Y.Value });
                }
                else
                {
                    Debug.Assert(point.Z != null, "point.Z != null");
                    coords = Coordinates(dim, lrs, geomTemplate,
                                         new[] { point.X.Value, point.Y.Value, point.Z.Value });
                }
                elemInfo = new double[] { 1, (int) SdoEType.Coordinate, 1 };
            }
            else
            {
                coords = Coordinates(dim, lrs, geomTemplate, ordinates);
            }

            switch ((SdoGTemplate)geomTemplate)
            {
                case SdoGTemplate.Coordinate:
                    return CreatePoint(factory, dim, lrs, elemInfo, 0, coords);

                case SdoGTemplate.Line:
                    return CreateLine(factory, dim, lrs, elemInfo, 0, coords);

                case SdoGTemplate.Polygon:
                    return CreatePolygon(factory, dim, lrs, elemInfo, 0, coords);

                case SdoGTemplate.MultiPoint:
                    return CreateMultiPoint(factory, dim, lrs, elemInfo, 0, coords);

                case SdoGTemplate.MultiLine:
                    return CreateMultiLine(factory, dim, lrs, elemInfo, 0, coords, -1);

                case SdoGTemplate.MultiPolygon:
                    return CreateMultiPolygon(factory, dim, lrs, elemInfo, 0, coords, -1);

                case SdoGTemplate.Collection:
                    return CreateCollection(factory, dim, lrs, elemInfo, 0, coords, -1);

                default:
                    return null;
            }
        }

        private static List<Coordinate> Coordinates(int dim, int lrs, int gtemplate, double[] ordinates)
        {
            if ((ordinates == null) || (ordinates.Length == 0))
            {
                return new List<Coordinate>();
            }

            //
            // POINT_TYPE Special Case
            //
            if ((dim == 2) && (lrs == 0) && (gtemplate == 01) && (ordinates.Length == 3))
            {
                var pt = new List<Coordinate>(1)
                {
                    new CoordinateZ(ordinates[0], ordinates[1], ordinates[2])
                };
                return pt;
            }

            int len = dim + lrs;

            if ((len == 0 && ordinates.Length != 0) || (len != 0 && ((ordinates.Length % len) != 0)))
            {
                throw new ArgumentException("Dimension D:" + dim + " and L:" +
                                         lrs + " denote Coordinates " + "of " + len +
                                         " ordinates. This cannot be resolved with" +
                                         "an ordinate array of length " + ordinates.Length);
            }

            int length = (len == 0 ? 0 : ordinates.Length / len);

            // we would have to ask for a dimension which represents all the requested
            // dimension and measures from a mask array in the future
            var pts = new List<Coordinate>(length);

            for (int i = 0; i < length; i++)
            {
                int offset = i * len;
                switch (len)
                {
                    case 2:
                        pts.Add(new CoordinateZ(ordinates[offset], ordinates[offset + 1], double.NaN));
                        break;
                    case 3:
                        pts.Add(new CoordinateZ(ordinates[offset], ordinates[offset + 1],
                                                ordinates[offset + 2]));
                        break;
                }

                //// in the future change this condition to include ignored dimensions from mask array
                //for (; j < actualDim && j < dim; j++)
                //{
                //    cs.setOrdinate(i, j, ordinates[i * len + j]);
                //    // may not always want to inc. j when we have a mask array
                //}
                ////// in the future change this condition to include ignored dimensions from mask array
                ////for (int d = j; j < actualDim && (j - d) < lrs; j++)
                ////{
                ////    cs.setOrdinate(i, j, ordinates[i * len + j]);
                ////    // may not always want to inc. j when we have a mask array
                ////}
            }
            return pts;
        }

        private GeometryCollection CreateCollection(GeometryFactory factory, int dim, int lrs, double[] elemInfo, int elemIndex,
                                                    List<Coordinate> coords, int numGeom)
        {

            int sOffset = StartingOffset(elemInfo, elemIndex);

            int length = coords.Count * dim;

            if (!(sOffset <= length))
                throw new ArgumentException("ELEM_INFO STARTING_OFFSET " + sOffset +
                                            " inconsistent with ORDINATES length " + coords.Count);

            int endTriplet = (numGeom != -1) ? elemIndex + numGeom : elemInfo.Length / 3 + 1;

            var list = new List<Geometry>();
            SdoEType etype;
            int interpretation;
            Geometry geom = null;

            bool cont = true;
            for (int i = elemIndex; cont && i < endTriplet; i++)
            {
                etype = EType(elemInfo, i);
                interpretation = Interpretation(elemInfo, i);

                switch (etype)
                {

                    case SdoEType.Unknown:
                        cont = false;
                        break;
                    case SdoEType.Coordinate:

                        if (interpretation == 1)
                        {
                            geom = CreatePoint(factory, dim, lrs, elemInfo, i, coords);
                        }
                        else if (interpretation > 1)
                        {
                            geom = CreateMultiPoint(factory, dim, lrs, elemInfo, i, coords);
                        }
                        else
                        {
                            throw new ArgumentException(
                                "ETYPE.POINT requires INTERPRETATION >= 1");
                        }

                        break;
                    case SdoEType.Line:
                        geom = CreateLine(factory, dim, lrs, elemInfo, i, coords);

                        break;

                    case SdoEType.Polygon:
                    case SdoEType.PolygonExterior:
                        geom = CreatePolygon(factory, dim, lrs, elemInfo, i, coords);
                        i += ((Polygon)geom).NumInteriorRings;

                        break;

                    case SdoEType.PolygonInterior:
                        throw new ArgumentException(
                            "ETYPE 2003 (Polygon Interior) no expected in a GeometryCollection" +
                         "(2003 is used to represent polygon holes, in a 1003 polygon exterior)");

                    default:
                        throw new ArgumentException("ETYPE " + etype +
                                                 " not representable as a JTS Geometry." +
                                                 "(Custom and Compound Straight and Curved Geometries not supported)");
                }
                if (cont)
                    list.Add(geom);
            }

            var geoms = factory.CreateGeometryCollection(list.ToArray());

            return geoms;
        }

        private MultiPolygon CreateMultiPolygon(GeometryFactory factory, int dim, int lrs, double[] elemInfo, int elemIndex,
                                                List<Coordinate> coords, int numGeom)
        {

            int sOffset = StartingOffset(elemInfo, elemIndex);
            var etype = EType(elemInfo, elemIndex);
            int interpretation = Interpretation(elemInfo, elemIndex);

            int length = coords.Count * dim;

            if (!(sOffset >= 1) || !(sOffset <= length))
                throw new ArgumentException("ELEM_INFO STARTING_OFFSET " + sOffset +
                                            " inconsistent with ORDINATES length " + coords.Count);

            if (etype != SdoEType.Polygon && etype != SdoEType.PolygonExterior)
                throw new ArgumentException("ETYPE " + etype + " inconsistent with expected POLYGON or POLYGON_EXTERIOR");

            if (interpretation != 1 && interpretation != 3)
            {
                return null;
            }

            int endTriplet = (numGeom != -1) ? elemIndex + numGeom : (elemInfo.Length / 3) + 1;

            var list = new List<Polygon>();
            bool cont = true;

            for (int i = elemIndex; cont && i < endTriplet && (etype = EType(elemInfo, i)) != SdoEType.Unknown; i++)
            {
                if ((etype == SdoEType.Polygon) || (etype == SdoEType.PolygonExterior))
                {
                    var poly = CreatePolygon(factory, dim, lrs, elemInfo, i, coords);
                    i += poly.NumInteriorRings; // skip interior rings
                    list.Add(poly);
                }
                else
                {
                    // not a Polygon - get out here
                    cont = false;
                }
            }

            var polys = factory.CreateMultiPolygon(list.ToArray());

            return polys;
        }

        private MultiLineString CreateMultiLine(GeometryFactory factory, int dim, int lrs, double[] elemInfo, int elemIndex,
                                                List<Coordinate> coords, int numGeom)
        {

            int sOffset = StartingOffset(elemInfo, elemIndex);
            var etype = EType(elemInfo, elemIndex);
            int interpretation = Interpretation(elemInfo, elemIndex);

            int length = coords.Count * dim;

            if (!(sOffset >= 1) || !(sOffset <= length))
                throw new ArgumentException("ELEM_INFO STARTING_OFFSET " + sOffset +
                                            " inconsistent with ORDINATES length " + coords.Count);
            if (!(etype == SdoEType.Line))
                throw new ArgumentException("ETYPE " + etype + " inconsistent with expected LINE");
            if (!(interpretation == 1))
            {
                // we cannot represent INTERPRETATION > 1
                return null;
            }

            int endTriplet = (numGeom != -1) ? (elemIndex + numGeom) : (elemInfo.Length / 3);

            var list = new List<LineString>();

            bool cont = true;
            for (int i = elemIndex; cont && i < endTriplet && (etype = EType(elemInfo, i)) != SdoEType.Unknown; i++)
            {
                if (etype == SdoEType.Line)
                {
                    list.Add(CreateLine(factory, dim, lrs, elemInfo, i, coords));
                }
                else
                {
                    // not a LineString - get out of here
                    cont = false;
                }
            }

            var lines = factory.CreateMultiLineString(list.ToArray());

            return lines;
        }

        private MultiPoint CreateMultiPoint(GeometryFactory factory, int dim, int lrs, double[] elemInfo, int elemIndex, List<Coordinate> coords)
        {
            int sOffset = StartingOffset(elemInfo, elemIndex);
            var etype = EType(elemInfo, elemIndex);
            int interpretation = Interpretation(elemInfo, elemIndex);

            if (!(sOffset >= 1) || !(sOffset <= coords.Count))
                throw new ArgumentException("ELEM_INFO STARTING_OFFSET " + sOffset +
                                            " inconsistent with ORDINATES length " + coords.Count);
            if (etype != SdoEType.Coordinate)
                throw new ArgumentException("ETYPE " + etype + " inconsistent with expected POINT");
            if (interpretation == 0)
            {
                return null;
            }

            int len = dim + lrs;

            int start = (sOffset - 1) / len;
            int end = start + interpretation;

            var points = factory.CreateMultiPointFromCoords(SubArray(coords, start, end));

            return points;
        }

        private Polygon CreatePolygon(GeometryFactory factory, int dim, int lrs, double[] elemInfo, int elemIndex, List<Coordinate> coords)
        {

            int sOffset = StartingOffset(elemInfo, elemIndex);
            var etype = EType(elemInfo, elemIndex);
            int interpretation = Interpretation(elemInfo, elemIndex);

            if (!(1 <= sOffset && sOffset <= (coords.Count * dim)))
            {
                throw new ArgumentException(
                    "ELEM_INFO STARTING_OFFSET " + sOffset +
                    "inconsistent with COORDINATES length " + (coords.Count * dim));
            }

            if (etype != SdoEType.Polygon && etype != SdoEType.PolygonExterior)
            {
                throw new ArgumentException("ETYPE " + etype + " inconsistent with expected POLYGON or POLYGON_EXTERIOR");
            }
            if (interpretation != 1 && interpretation != 3)
            {
                return null;
            }

            var exteriorRing = CreateLinearRing(factory, dim, lrs, elemInfo, elemIndex, coords);

            var rings = new List<LinearRing>();

            bool cont = true;
            for (int i = elemIndex + 1; cont && (etype = EType(elemInfo, i)) != SdoEType.Unknown; i++)
            {
                if (etype == SdoEType.PolygonInterior)
                {
                    rings.Add(CreateLinearRing(factory, dim, lrs, elemInfo, i, coords));
                }
                else if (etype == SdoEType.Polygon)
                {
                    // need to test Clockwiseness of Ring to see if it is
                    // interior or not - (use POLYGON_INTERIOR to avoid pain)

                    var ring = CreateLinearRing(factory, dim, lrs, elemInfo, i, coords);

                    if (Algorithm.Orientation.IsCCW(ring.CoordinateSequence))
                    {
                        // it is an Interior Hole
                        rings.Add(ring);
                    }
                    else
                    {
                        // it is the next Polygon! - get out of here
                        cont = false;
                    }
                }
                else
                {
                    // not a LinearRing - get out of here
                    cont = false;
                }
            }

            var poly = factory.CreatePolygon(exteriorRing, rings.ToArray());

            return poly;
        }


        private LinearRing CreateLinearRing(GeometryFactory factory, int dim, int lrs, double[] elemInfo, int elemIndex, List<Coordinate> coords)
        {

            int sOffset = StartingOffset(elemInfo, elemIndex);
            var etype = EType(elemInfo, elemIndex);
            int interpretation = Interpretation(elemInfo, elemIndex);
            int length = coords.Count * dim;

            if (!(sOffset <= length))
                throw new ArgumentException("ELEM_INFO STARTING_OFFSET " + sOffset +
                                            " inconsistent with ORDINATES length " + coords.Count);
            if (etype != SdoEType.Polygon && etype != SdoEType.PolygonExterior &&
                etype != SdoEType.PolygonInterior)
            {
                throw new ArgumentException("ETYPE " + etype +
                                            " inconsistent with expected POLYGON, POLYGON_EXTERIOR or POLYGON_INTERIOR");
            }
            if (interpretation != 1 && interpretation != 3)
            {
                return null;
            }
            LinearRing ring;

            int len = (dim + lrs);
            int start = (sOffset - 1) / len;
            int eOffset = StartingOffset(elemInfo, elemIndex + 1); // -1 for end
            int end = (eOffset != -1) ? ((eOffset - 1) / len) : coords.Count;

            if (interpretation == 1)
            {
                ring = new LinearRing(ToPointArray(SubList(coords, start, end)));
            }
            else
            {
                // interpretation == 3
                // rectangle does not maintain measures
                var pts = new List<Coordinate>(5);
                var ptssrc = SubList(coords, start, end);
                var min = ptssrc[0];
                var max = ptssrc[1];
                pts.AddRange(new[]
                                 {
                                     min, new Coordinate(max.X, min.Y), max, new Coordinate(min.X, max.Y) , min
                                 });

                ring = factory.CreateLinearRing(pts.ToArray());
            }

            return ring;
        }


        private LineString CreateLine(GeometryFactory factory, int dim, int lrs, double[] elemInfo, int elemIndex, List<Coordinate> coords)
        {

            int sOffset = StartingOffset(elemInfo, elemIndex);
            var etype = EType(elemInfo, elemIndex);
            int interpretation = Interpretation(elemInfo, elemIndex);

            if (etype != SdoEType.Line)
                return null;

            if (interpretation != 1)
            {
                throw new ArgumentException("ELEM_INFO INTERPRETAION " +
                                         interpretation + " not supported" +
                                         "by JTS LineString.  Straight edges" +
                                         "( ELEM_INFO INTERPRETAION 1) is supported");
            }

            int
        len = (dim + lrs);
            int start = (sOffset - 1) / len;
            int eOffset = StartingOffset(elemInfo, elemIndex + 1); // -1 for end
            int end = (eOffset != -1) ? ((eOffset - 1) / len) : coords.Count;


            var line = factory.CreateLineString(ToPointArray(SubList(coords, start, end)));

            return line;
        }



        private static Coordinate[] ToPointArray(ICollection<Coordinate> input)
        {
            var pts = new List<Coordinate>(input.Count);
            foreach (var point in input)
                pts.Add(new CoordinateZ(point.X, point.Y, point.Z));

            return pts.ToArray();
        }

        private Point CreatePoint(GeometryFactory factory, int dim, int lrs, double[] elemInfo, int elemIndex, List<Coordinate> coords)
        {
            int sOffset = StartingOffset(elemInfo, elemIndex);
            var etype = EType(elemInfo, elemIndex);
            int interpretation = Interpretation(elemInfo, elemIndex);

            if (!(sOffset >= 1) || !(sOffset <= coords.Count))
                throw new ArgumentException("ELEM_INFO STARTING_OFFSET " + sOffset +
                                            " inconsistent with ORDINATES length " + coords.Count);
            if (etype != SdoEType.Coordinate)
                throw new ArgumentException("ETYPE " + etype + " inconsistent with expected POINT");
            if (interpretation != 1)
            {
                return null;
            }

            int len = (dim + lrs);
            int start = (sOffset - 1) / len;
            int eOffset = StartingOffset(elemInfo, elemIndex + 1); // -1 for end

            Coordinate point;
            if ((sOffset == 1) && (eOffset == -1))
            {
                // Use all Coordinates
                point = coords[0];
            }
            else
            {
                int end = (eOffset != -1) ? ((eOffset - 1) / len) : coords.Count;
                point = SubList(coords, start, end)[0];
            }

            return factory.CreatePoint(point);
        }



        private static List<Coordinate> SubList(List<Coordinate> coords, int start, int end)
        {
            if ((start == 0) && (end == coords.Count))
            {
                return coords;
            }

            return coords.GetRange(start, (end - start));
        }

        private static Coordinate[] SubArray(List<Coordinate> coords, int start, int end)
        {
            return coords.GetRange(start, end - start).ToArray();
        }


        private static SdoEType EType(double[] elemInfo, int tripletIndex)
        {
            if (((tripletIndex * 3) + 1) >= elemInfo.Length)
            {
                return SdoEType.Unknown;
            }

            return (SdoEType)elemInfo[(tripletIndex * 3) + 1];
        }

        private static int Interpretation(double[] elemInfo, int tripletIndex)
        {
            if (((tripletIndex * 3) + 2) >= elemInfo.Length)
            {
                return -1;
            }

            return (int) elemInfo[(tripletIndex * 3) + 2];
        }

        private static int StartingOffset(double[] elemInfo, int tripletIndex)
        {
            if (((tripletIndex * 3) + 0) >= elemInfo.Length)
            {
                return -1;
            }

            return (int) elemInfo[(tripletIndex * 3) + 0];
        }
    }
}
