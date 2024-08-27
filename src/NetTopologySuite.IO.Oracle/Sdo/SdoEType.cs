namespace NetTopologySuite.IO.Sdo
{
#pragma warning disable 1591
    internal enum SdoEType
    {
        Unknown = -1,

        Coordinate = 1,
        Line = 2,
        Polygon = 3,

        PolygonExterior = 1003,
        PolygonInterior = 2003
    }
#pragma warning restore 1591
}
