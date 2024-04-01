using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
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

        public void CalculateValue(Expression exp)
        {
            value = exp.CalculateValue(coords);
        }
    }
}
