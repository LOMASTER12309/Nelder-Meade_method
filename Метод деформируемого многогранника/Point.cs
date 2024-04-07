using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;

using FunctionParser;

namespace Метод_деформируемого_многогранника
{
    class ValueIncreasingComparer : IComparer<Point>
    {
        public int Compare(Point ind1, Point ind2)
        {
            return Math.Sign(ind1.value - ind2.value);
        }
    }
    class ValueDecreasingComparer : IComparer<Point>
    {
        public int Compare(Point ind1, Point ind2)
        {
            return -Math.Sign(ind1.value - ind2.value);
        }
    }
    class Point
    {
        public int dimension = 0;
        public double[] coords;
        public double value = double.NaN;

        public Point(double[] arr)
        {
            coords = (double[])arr.Clone();
            dimension = coords.Length;
        }
        public Point(double[] arr, Expression exp) : this(arr)
        {
            if (exp != null) value = exp.CalculateValue(coords);
        }
        public void PrintPoint(int precision = 3, string name = null)
        {
            Console.Write($"{name} = ");
            precision = Math.Abs(precision);
            Console.Write("(");
            for (int i = 0; i < dimension - 1; i++) Console.Write($"{Math.Round(coords[i], precision)}, ");
            Console.Write($"{Math.Round(coords[dimension - 1], precision)})");
        }
        public void PrintPointWithValue(string name = null, int precision = 3)
        {
            PrintPoint(precision, name);
            precision = Math.Abs(precision);
            if (name == null)
            {
                name = "";
                for (int i = 0; i < dimension - 1; i++) name += $"{Math.Round(coords[i], precision)}, ";
                name += $"{Math.Round(coords[dimension - 1], precision)}";
            }
            Console.Write($", F({name}) = {value}");
        }
        public static Point operator -(Point point1, Point point2)
        {
            if (point1.dimension != point2.dimension)
                throw new Exception("Не совпадают размерности точек");
            double[] newCoords = new double[point1.dimension];
            for (int i = 0; i < point1.dimension; i++)
                newCoords[i] = point1.coords[i] - point2.coords[i];
            return new Point(newCoords);
        }

        public static Point operator +(Point point1, Point point2)
        {
            if (point1.dimension != point2.dimension)
                throw new Exception("Не совпадают размерности точек");
            double[] newCoords = new double[point1.dimension];
            for (int i = 0; i < point1.dimension; i++)
                newCoords[i] = point1.coords[i] + point2.coords[i];
            return new Point(newCoords);
        }
        public static bool LinearIndependence(Point point1, Point point2)
        {
            if (point1.dimension != point2.dimension)
                throw new Exception("Не совпадают размерности точек");
            double[] coords1 = point1.coords;
            double[] coords2 = point2.coords;
            double ratio = 0;
            for (int i = 0; i < point1.dimension; i++)
                if ((coords1[i] != 0) && (coords2[i] != 0))
                {
                    ratio = coords1[i] / coords2[i];
                    break;
                }
            for (int i = 0; i < point1.dimension; i++)
            {
                if ((coords1[i] == 0) && (coords2[i] == 0)) continue;
                else if ((coords1[i] * coords2[i] == 0) && (coords2[i] + coords2[i] != 0)) return true;
                else if (coords1[i] / coords2[i] != ratio) return true;
            }
            return false;
        }
        public void CalculateValue(Expression exp)
        {
            value = exp.CalculateValue(coords);
        }
        public static void CalculateValue(Point[] points, Expression exp)
        {
            foreach (var point in points) point.CalculateValue(exp);
        }
        public static Point Clone(Point point)
        {
            Point newPoint = new Point(point.coords);
            newPoint.value = point.value;
            return newPoint;
        }
        public static Point[] Clone(Point[] points)
        {
            Point[] newPoints = new Point[points.Length];
            for (int i = 0; i < points.Length; i++)
                newPoints[i] = Clone(points[i]);
            return newPoints;
        }
    }
}
