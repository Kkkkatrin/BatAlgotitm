using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace BatAlgorithmC
{
    class Program
    {
        static void Main(string[] args)
        {
                
                Mybat x = new Mybat(200, 100, 0.5, 0.7, 4, Mybat.Result);
                Console.WriteLine("Hit any key to continue.");
                Console.ReadLine();

        }


    }

        public class Mybat
        {

            public int _numberOfBats, _generations, Qmin, Qmax, N_iter, _dimension;
            public double _volume, _pulseRate, min, max, fnew, fmin;
            public double[][] _lowerBound, _upperBound, _velocity, _solution, S, velocity;
            public double[] _fitness, _tempSolution, _bestSolution, Q;
            public Random random;
           
            public static void initJagged(double[][] array, int n, int d)
            {
                for (int i = 0; i < n; i++) array[i] = new double[d];
            }

            public Mybat(
                int bats,
                int generations,
                double loud,
                double pulse,
                int dimension,
                Func<double[], int, double> function
                )
            {
                //инициализация переменных
                _numberOfBats = bats;
                _generations = generations;
                _volume = loud;
                _pulseRate = pulse;
                _dimension = dimension;
                
               
                Random random = new Random();

                min = -15;
                max = 15;

                fmin = 0;
                //спуск за пределы
                _lowerBound = new double[1][];
                _upperBound = new double[1][];

                Q = new double[_numberOfBats]; // частота
                _velocity = new double[_numberOfBats][]; //скорость
                velocity = new double[_numberOfBats][];
               
                initJagged(_velocity, _numberOfBats, _dimension);
                initJagged(_lowerBound, 1, _dimension);
                initJagged(_upperBound, 1, _dimension);


                //Иннициализация массивов решений 
                _solution = new double[_numberOfBats][];
                S = new double[_numberOfBats][];
                _fitness = new double[_numberOfBats]; // временный контейнер
                _bestSolution = new double[_dimension];
                _tempSolution = new double[_dimension]; //Временный держатель строки в массиве

                initJagged(_solution, _numberOfBats, _dimension);
                initJagged(S, _numberOfBats, _dimension);



                //for (int i = 0; i < _numberOfBats; i++)
                //{
                //    // Для уменьшения кол-ва кода добавлен массив Q[]array с '0'
                //    Q[i] = 0;
                
                //    for (int x = 0; x < _dimension; x++)
                //    {
                //        // Для уменьшения кол-ва кода добавлен массив initialize _velocity[][] array with '0' as element
                //        _velocity[i][x] = 0;

                //        //поиск случайных двойных значений от lowerBound до upperBound
                //        _solution[i][x] = (random.NextDouble() * (max - min)) + min;
                //        _tempSolution[x] = _solution[i][x];
                //    }
                //    _fitness[i] = function(_tempSolution, _dimension);

                //    //Инициализация лучшего из минимумов
                //    if (i == 0 || fmin > _fitness[i])
                //    {
                //        fmin = _fitness[i];

                //        for (int x = 0; x < _dimension; x++)
                //        {
                //            _bestSolution[x] = _solution[i][x];
                //        }
                //    }
                //    Console.WriteLine("fitness[" + i + "]" + _fitness[i]); //test
                //}
                //Console.WriteLine("fmin = " + fmin); //test

                // при необходимости изменить для максимальной эффективности
                Qmin = 0;
                Qmax = 5;
                N_iter = 100; //количество оценок функции

               
               
                // собственно та самая летучая мышь
                for (int loop = 0; loop < N_iter; loop++)
                {
                    //перебрать все биты / решения
                    for (int nextBat = 0; nextBat < _numberOfBats; nextBat++)
                    {
                        Q[nextBat] = Qmin + ((Qmin - Qmax) * random.NextDouble());
                        

                        // цикл для ускорения процесса
                        for (int vel = 0; vel < _dimension; vel++)
                        {
                            _velocity[nextBat][vel] = _velocity[nextBat][vel] + ((_solution[nextBat][vel] - _bestSolution[vel]) * Q[nextBat]);
                        }

                        //новые решения
                        for (int nextDimension = 0; nextDimension < _dimension; nextDimension++)
                        {
                            S[nextBat][nextDimension] = _solution[nextBat][nextDimension] + _velocity[nextBat][nextDimension];
                        }

                        // частота пульса
                        if (random.NextDouble() > _pulseRate)
                        {
                            for (int nextDimension = 0; nextDimension < _dimension; nextDimension++)
                            {
                                S[nextBat][nextDimension] = _bestSolution[nextDimension] + (0.001 * random.NextGaussian());
                            }
                        }

                        //помещение текущей строки во временный массив
                        for (int nextDimension = 0; nextDimension < _dimension; nextDimension++)
                        {
                            _tempSolution[nextDimension] = S[nextBat][nextDimension];
                        }
                        fnew = function(_tempSolution, _dimension);

                        // обновить, если решение улучшено, и не слишком громко
                        if ((fnew <= _fitness[nextBat]) && (random.NextDouble() < _volume))
                        {
                            for (int x = 0; x < _dimension; x++)
                            {
                                _solution[nextBat][x] = S[nextBat][x];
                                _fitness[nextBat] = fnew;
                            }
                        }

                        //обновить текущее лучшее решение
                        if (fnew <= fmin)
                        {
                            for (int nextDimension = 0; nextDimension < _dimension; nextDimension++)
                            {
                                _bestSolution[nextDimension] = S[nextBat][nextDimension];
                                fmin = fnew;
                            }
                        }
                    }
                }

                Console.WriteLine(" ");
                for (int i = 0; i < _numberOfBats; i++)
                {
                    Console.WriteLine("fitness[" + i + "]" + _fitness[i]);
                }

                for (int nextDimension = 0; nextDimension < _dimension; nextDimension++)
                {
                    Console.WriteLine("best[" + nextDimension + "]" + _bestSolution[nextDimension]);
                }
                Console.WriteLine("Fmin = " + fmin);
            
                for (int i = 0; i < _numberOfBats; i++)
                {
                    // Для уменьшения кол-ва кода добавлен массив Q[]array с '0'
                    Q[i] = 0;

                    for (int x = 0; x < _dimension; x++)
                    {
                        // Для уменьшения кол-ва кода добавлен массив initialize _velocity[][] array with '0' as element
                        _velocity[i][x] = 0;

                        //поиск случайных двойных значений от lowerBound до upperBound
                        _solution[i][x] = (random.NextDouble() * (max - min)) + min;
                        _tempSolution[x] = _solution[i][x];
                    }
                    _fitness[i] = function(_tempSolution, _dimension);

                    //Инициализация лучшего из минимумов
                    if (i == 0 || fmin > _fitness[i])
                    {
                        fmin = _fitness[i];

                        for (int x = 0; x < _dimension; x++)
                        {
                            _bestSolution[x] = _solution[i][x];
                        }
                    }
                    Console.WriteLine("fitness[" + i + "]" + _fitness[i]); //test
                }
                Console.WriteLine("fmin = " + fmin); //test
            }

            public void set_bounds(int x, double L, double U)
            {
                //double temp_Lb[x];
                //double temp_Ub[x];
                for (int i = 0; i < x; i++)
                {
                    _lowerBound[0][i] = L;
                    _upperBound[0][i] = U;
                }
            }
            public static double Result(double[] value, int d) 
            {
                //тестовая функция Стыбинского-Танга
                double result = 0;
                
                for (int i = 0; i < d; i++)
                {
                    result = (Math.Pow(value[i], 4) - 16 * Math.Pow(value[i], 2) + 5 * value[i]) / 2;
                }
                return result;
            }


            public static double sphere(double[] value, int d)
            {
                // сферическая функция, где fmin находится в 0
                double result = 0;

                for (int i = 0; i < d; i++)
                {
                    result += (value[i] * value[i]);
                }
                return result;
            }

            public static double squareRoot(double[] value, int d)
            {
                // нахождение квадратного корня для 0.3
                double result = 0;

                for (int i = 0; i < d; i++)
                {
                    result += Math.Abs(.3 - (value[i] * value[i]));
                }
                return result;
            }

            public static double RootOfXYEquations(double[] value, int d)
            {
                // решение для x и y  xy = 6  and x+y = 5
                double result = 0;

                result += Math.Abs(5 - (value[0] + value[1]));
                result += Math.Abs(6 - (value[0] * value[1]));


                return result;
            }
        }
        

        static class MathExtensiionns
        {
            public static double NextGaussian(this Random rand)
            {
                double u1 = rand.NextDouble(); //это равномерные (0,1) случайные двойники
                double u2 = rand.NextDouble();
                double mean = 0, stdDev = 1;
                double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                             Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
                double randNormal =
                             mean + stdDev * randStdNormal; //random normal(mean,stdDev^2)
                return randNormal;
            }
        }
    }
