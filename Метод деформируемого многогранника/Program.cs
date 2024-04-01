using System.Drawing;
using FunctionParser;

namespace Метод_деформируемого_многогранника
{
    internal class Program
    {
        static void Main(string[] args)
        {
            int n = 0; //размерность
            Expression exp = new Expression("0", null, null);
            double alpha = 0; //коэффициент отражения
            double beta = 0;  //коэффициент растяжения
            double gamma = 0; //коэффициент сжатия
            double l = 0; //начальное отклонение
            double epsilon = 0; //эпсилон
            bool update = false;
            int FrequencyOfUpdates = -1; //период обновления симплекса
            int precision = -1; //точность округление
            Console.Write("Количество переменных (размерность) n = ");
            while (n <= 0) n = protect_int_input();
            var arg = new string[n];
            for (int i = 0; i < n; i++) arg[i] = $"x{i + 1}";
            Console.Write("Функция от переменных ");
            for (int i = 0; i < n - 1; i++) Console.Write($"{arg[i]}, ");
            Console.Write($"{arg[n - 1]}\n f(X) = ");
            string expression;
            bool property = false;
            do
            {
                expression = Console.ReadLine();
                if (!(property = Expression.IsExpression(expression, arg)))
                    Console.Write("Некорректное выражение! Попробуйте ещё: f(x) = ");
            }
            while (!property);
            exp = new Expression(expression, arg, null);
            Console.Write("Коэффициент отражения alpha = ");
            while (alpha <= 0) alpha = protect_double_input();
            Console.Write("Коэффициент растяжения beta = ");
            while (beta <= 1) beta = protect_double_input();
            Console.Write("Коэффициент сжатия gamma = ");
            while ((gamma >= 1) || (gamma <= 0)) gamma = protect_double_input();
            Console.Write("Начальное отклонение l = ");
            while (l <= 0) l = protect_double_input();
            Console.Write("Для условия останова epsilon = ");
            while (epsilon <= 0) epsilon = protect_double_input();
            string c = " ";
            Console.Write("Производить периодическое обновление симплекса? (1-да, 0-нет): ");
            while ((c != "1") && (c != "0")) c = Console.ReadLine();
            if (c == "1")
            {
                update = true;
                Console.Write("Период обновления симлпекса FrequencyOfUpdates = ");
                while (FrequencyOfUpdates < 2) FrequencyOfUpdates = protect_int_input();
            }
            Console.WriteLine("Округление (число знаков после запятой) = ");
            while (precision < 0) precision = protect_int_input();
            Console.Write("Начальная точка X0 = ");
            string[] x0;
            double[] StatrPoint;
            Point[] simplex = new Point[n + 1];
            bool correct = false;
            do
            {
                x0 = Console.ReadLine().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
                if (x0.Length < n)
                {
                    Console.Write("Неверный формат! X0 = ");
                    continue;
                }
                StatrPoint = new double[n];
                for (int i = 0; i < n; i++)
                    if (!double.TryParse(x0[i], out StatrPoint[i]))
                    {
                        Console.Write("Неверный формат! X0 = ");
                        continue;
                    }
                simplex[0] = new Point(StatrPoint, exp);
                correct = true;
            }
            while (!correct);
            GenerateSimplex(simplex, exp, l);
            Array.Sort(simplex, new ValueIncreasingComparer());
            Console.WriteLine("Начальный сиплекс: ");
            PrintSimplex(simplex, 1, precision);
            int k = 0;
            Point centrePoint;
            Point refPoint;
            Point stretchPoint;
            Point compressPoint;
            double dispersion = Dispersion(simplex);
            Console.WriteLine($"Дисперсия = {dispersion}");
            Console.WriteLine("Начало работы алгоритма: ");

            while (dispersion > epsilon)
            {
                k++;
                Console.WriteLine($"Цикл {k}:");
                if ((k % FrequencyOfUpdates == 0) && (k > 1) && (update))
                {
                    Console.WriteLine("Обновление симплекса!");
                    GenerateSimplex(simplex, exp, l);
                    Array.Sort(simplex, new ValueIncreasingComparer());
                    Console.WriteLine("Новый симплекс:");
                    PrintSimplex(simplex, 1, precision);
                    Console.WriteLine($"Дисперсия = {dispersion = Dispersion(simplex)}");
                    Console.WriteLine();
                    continue;
                }
                centrePoint = Centre(simplex);
                refPoint = Reflection(simplex, centrePoint, alpha, exp);
                refPoint.PrintPointWithValue("U(k)", precision);
                Console.WriteLine();
                if ((simplex[0].value <= refPoint.value) && (refPoint.value <= simplex[n - 1].value)) //случай 1
                {
                    simplex[n] = refPoint;
                    Console.WriteLine($"Заменяем Xn на U(k)");
                }
                else if (refPoint.value < simplex[0].value) //случай 2
                {
                    stretchPoint = Stretching(centrePoint, refPoint, beta, exp);
                    stretchPoint.PrintPointWithValue("V(k)", precision);
                    Console.WriteLine();
                    if (stretchPoint.value < refPoint.value)
                    {
                        simplex[n] = stretchPoint;
                        Console.WriteLine($"Заменяем Xn на V(k)");
                    }
                    else
                    {
                        simplex[n] = refPoint;
                        Console.WriteLine($"Заменяем Xn на U(k)");
                    }
                }
                else if (refPoint.value > simplex[n - 1].value) //случай 3
                {
                    compressPoint = Compression(simplex, centrePoint, refPoint, gamma, exp);
                    compressPoint.PrintPointWithValue("W(k)", precision);
                    Console.WriteLine();
                    if (compressPoint.value < Math.Min(simplex[n].value, refPoint.value))
                    {
                        simplex[n] = compressPoint;
                        Console.WriteLine($"Заменяем Xn на W(k)");
                    }
                    else
                    {
                        Console.WriteLine("Производим глобальное сжатие симплекса");
                        GlobalComprerssion(simplex, exp);
                    }
                }
                Array.Sort(simplex, new ValueIncreasingComparer());
                Console.WriteLine("Новый симплекс:");
                PrintSimplex(simplex, 1, precision);
                Console.WriteLine($"Дисперсия = {dispersion = Dispersion(simplex)}");
                Console.WriteLine();
            }
            Console.Write("Лучший результат: ");
            simplex[0].PrintPointWithValue("X0", precision);
        }
        static double Dispersion(Point[] simplex)
        {
            int n = simplex[0].dimension;
            double disp = 0;
            for (int j = 1; j < n + 1; j++)
                disp += Math.Pow(simplex[j].value - simplex[0].value, 2);
            disp /= n;
            return Math.Sqrt(disp);
        }
        static Point Centre(Point[] simplex)
        {
            int n = simplex[0].dimension;
            double[] point = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    point[i] += simplex[j].coords[i];
                point[i] /= n;
            }
            return new Point(point);
        }
        static Point Reflection(Point[] simplex, Point centre, double alpha, Expression exp)
        {
            int n = simplex[0].dimension;
            double[] point = new double[n];
            for (int i = 0; i < n; i++)
                point[i] = centre.coords[i] + alpha * (centre.coords[i] - simplex[n].coords[i]);
            return new Point(point, exp);
        }
        static Point Stretching(Point centre, Point refPoint, double beta, Expression exp)
        {
            int n = centre.dimension;
            double[] point = new double[n];
            for (int i = 0; i < n; i++)
                point[i] = centre.coords[i] + beta * (refPoint.coords[i] - centre.coords[i]);
            return new Point(point, exp);
        }

        static Point Compression(Point[] simplex, Point centre, Point refPoint, double gamma, Expression exp)
        {
            int n = simplex[0].dimension;
            double[] point = new double[n];
            if (simplex[n].value <= refPoint.value)
                for (int i = 0; i < n; i++)
                    point[i] = centre.coords[i] + gamma * (simplex[n].coords[i] - centre.coords[i]);
            else
                for (int i = 0; i < n; i++)
                    point[i] = centre.coords[i] + gamma * (refPoint.coords[i] - centre.coords[i]);
            return new Point(point, exp);
        }

        static void GlobalComprerssion(Point[] simplex, Expression exp)
        {
            int n = simplex[0].dimension;
            double[] newPoint;
            for (int i = 1; i < n + 1; i++)
            {
                newPoint = new double[n];
                for (int j = 0; j < n; j++)
                    newPoint[j] = 0.5 * (simplex[i].coords[j] + simplex[0].coords[j]);
                simplex[i] = new Point(newPoint, exp);
            }
        }
        static void GenerateSimplex(Point[] simplex, Expression exp, double l)
        {
            int n = simplex[0].dimension;
            double[] coords;
            for (int i = 1; i < n + 1; i++)
            {
                coords = (double[])simplex[0].coords.Clone();
                coords[i - 1] += l;
                simplex[i] = new Point(coords, exp);
            }
        }

        static void PrintSimplex(Point[] simplex, int NumInLine = 1, int precision = 3)
        {
            int n = simplex[0].dimension;
            for (int i = 0; i < n; i++)
            {
                simplex[i].PrintPointWithValue($"X{i}");
                Console.Write("; ");
                if ((i + 1) % NumInLine == 0) Console.WriteLine();
            }
            simplex[n].PrintPointWithValue($"X{n}");
            Console.WriteLine();
        }
        static double protect_double_input()
        {
            double num = 0;
            while (true)
            {
                string text = Console.ReadLine();
                if (double.TryParse(text, out num)) break;
                Console.Write("Неверный формат! Попробуйте ещё: ");
            }
            return num;
        }
        static int protect_int_input()
        {
            int num = 0;
            while (true)
            {
                string text = Console.ReadLine();
                if (int.TryParse(text, out num)) break;
                Console.Write("Неверный формат! Попробуйте ещё: ");
            }
            return num;
        }
    }
}
