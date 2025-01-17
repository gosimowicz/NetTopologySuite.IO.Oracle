using Oracle.ManagedDataAccess.Client;
using System.Configuration;
using System.Data;
using NUnit.Framework;
using NetTopologySuite.IO.Sdo;
using System;

namespace NetTopologySuite.IO.Oracle.Connection.Test
{
    public static class OracleHelper
    {
        /// <summary>
        /// Opens a connection to the test database
        /// </summary>
        /// <returns></returns>
        public static OracleConnection OpenConnection(string connectionString = null)
        {
            connectionString ??= ConfigurationManager.AppSettings.Get("TestDBConnectionString");
            try
            {
                if (string.IsNullOrWhiteSpace(connectionString))
                {
                    Assert.Ignore("Connection string is empty!");
                    return null;
                }
                //TestContext.Error.WriteLine("Trying to connect with '{0}'", connectionString);
                var conn = new OracleConnection(connectionString);
                conn.Open();
                //TestContext.Error.WriteLine("Connection successful!");
                //TestContext.Error.WriteLine("Connected to '{0}' on '{1}'.", conn.DatabaseName, conn.DatabaseEditionName);
                return conn;
            }
            catch (Exception ex)
            {
                TestContext.Error.WriteLine(ex.Message);
                TestContext.Error.WriteLine(ex.StackTrace);
                throw new IgnoreException("Connection to Oracle database server failed", ex);
            }
        }

        /// <summary>
        /// Drops (if it exists) and recreates the given table.
        /// </summary>
        /// <param name="connection"></param>
        /// <param name="tableName"> </param>
        /// <returns>The name of the table as returned from sys.all_tables</returns>
        public static string CreateGeometryTable(OracleConnection connection, string tableName)
        {
            DropGeometryTable(connection, tableName);

            string queryString = $"CREATE TABLE {tableName} (data MSYS.SDO_GEOMETRY)";
            using var command = new OracleCommand(queryString, connection);
            command.ExecuteNonQuery();

            string queryString2 = $"SELECT TABLE_NAME FROM sys.all_tables WHERE TABLE_NAME = '{tableName}'";
            using var command2 = new OracleCommand(queryString2, connection);
            return (string)command2.ExecuteScalar();
        }

        /// <summary>
        /// Drops (if it exists) the given table.
        /// </summary>
        /// <param name="connection"></param>
        /// <param name="tableName"></param>
        public static void DropGeometryTable(OracleConnection connection, string tableName)
        {
            // Remove 'GEO_DATA' table
            string queryString = $@"BEGIN
                EXECUTE IMMEDIATE 'DROP TABLE {tableName}';
                EXCEPTION
                    WHEN OTHERS THEN
                        IF SQLCODE != -942 THEN
                            RAISE;
                        END IF;
                END;";
            using var command = new OracleCommand(queryString, connection);
            command.ExecuteNonQuery();
        }


        /// <summary>
        /// Assumption GEO_DATA table exists.
        /// Write a new created Geometry object to the database.
        /// </summary>
        public static Geometries.Geometry WriteGeometryToTable(OracleConnection connection, string wkt, string testTableName)
        {
            var geom = ConvertWKTToGeometry(wkt);
            var udt = ConvertGeomToOracleUDT(geom);

            // Open the connection
            bool wasClosed = connection.State == ConnectionState.Closed;
            if (wasClosed) connection.Open();

            // Drop & Create Geometry table.
            CreateGeometryTable(connection, testTableName);

            string queryString = $"INSERT INTO {testTableName} (data) VALUES (:geo)";

            using var command = new OracleCommand(queryString, connection);
            var geometryParam = new OracleParameter()
            {
                ParameterName = "geo",
                DbType = DbType.Object,
                Value = udt,
                Direction = ParameterDirection.Input,
                UdtTypeName = "MDSYS.SDO_GEOMETRY"
            };
            command.Parameters.Add(geometryParam);
            command.ExecuteNonQuery();

            // Close connection
            if (wasClosed) connection.Close();

            return geom;
        }

        private static Sdo.SdoGeometry ConvertGeomToOracleUDT(Geometries.Geometry geom)
        {
            // Write geometry object into UDT object.
            var oracleWriter = new OracleGeometryWriter();
            var udt = oracleWriter.Write(geom);
            return udt;
        }

        private static Geometries.Geometry ConvertWKTToGeometry(string wkt)
        {
            // Read WKT into geometry object.
            var wr = new WKTReader { IsOldNtsCoordinateSyntaxAllowed = false };
            string correctCCW = wkt;
            var geom = wr.Read(correctCCW);
            return geom;
        }

        /// <summary>
        /// Read a newly created Geometry object from the database.
        /// Assumption GEO_DATA table exists.
        /// </summary>
        public static Geometries.Geometry ReadGeometryFromTable(OracleConnection connection, string testTableName)
        {
            // Open connection
            bool wasClosed = connection.State == ConnectionState.Closed;
            if (wasClosed) connection.Open();

            // Write query string & command
            string queryString = $"SELECT * FROM {testTableName}";
            using var command = new OracleCommand(queryString, connection);

            var geometryParam = new OracleParameter();
            command.Parameters.Add(geometryParam);
            var res = (SdoGeometry) command.ExecuteScalar();

            var oracleReader = new OracleGeometryReader(NtsGeometryServices.Instance);
            var geom2 = oracleReader.Read(res);

            // Close connection
            if (wasClosed) connection.Close();

            return geom2;
        }


    }
}
