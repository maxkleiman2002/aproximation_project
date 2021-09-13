#include <iostream>
#include <math.h>
#include <Windows.h>
#include <vector>
#include <fstream>
#include <cstdlib>

using namespace std;


class Approximation {

private:
    int m_n;
    int m_k;
    int i;
    int n;


public:


    void print_array() {
        std::vector<double> x;
        std::vector<double> y;

        ifstream inf("info_in.txt");
        ifstream infn("info_in_y.txt");
        ofstream outf("info_out.txt");
        if (!inf) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else if (!infn) {
            cerr << "Не удалось считать файл # 2! " << endl;
            exit(1);
        }
        else {
            int n;

            inf >> n;
            outf << n << endl;
            for (int i = 0; i < n; i++) {
                double k;
                inf >> k;
                x.push_back(k);
                outf << "X[" << i << "] = ";
                outf << "\t" << x[i];
                outf << endl;
                outf << endl;
            }

            for (int j = 0; j < n; j++) {
                double t;
                infn >> t;
                y.push_back(t);
                outf << "Y[" << j << "] = ";
                outf << "\t" << y[j];
                outf << endl;
                outf << endl;

            }


            if (!outf) {
                cerr << "Не удалось записать в  файл! " << endl;
                exit(1);
            }
            else
            {
                cout << "Массивы чисел X и Y записаны в файл." << endl;

                cout << " X = ";
                for (int i = 0; i < n; i++) {
                    cout << "\t" << x[i];
                }
                cout << endl;

                cout << " Y = ";
                //outf << "Y = ";
                for (int j = 0; j < n; j++) {
                    cout << "\t" << y[j];
                }
                cout << endl;

            }

        }
    }
    void Aprox_linear_func() {
        vector<double> x;
        vector<double> y;

        double sx = 0, sy = 0, sxy = 0, sxx = 0;
        ifstream inf("info_in.txt");
        ifstream infn("info_in_y.txt");
        ofstream outf("info_out.txt");
        if (!inf) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else {
            int n;
            inf >> n;
            double x0 = 1.2;
            for (int i = 0; i < n; i++) {
                double k;
                inf >> k;
                x.push_back(k);
                double t;
                infn >> t;
                y.push_back(t);
                sx += x[i];
                sy += y[i];
                sxx += x[i] * x[i];
                sxy += x[i] * y[i];

            }
            double K = (sx * sy - n * sxy) / (sx * sx - n * sxx);
            double B = (sy - K * sx) / n;
            double R = K * x0 + B;
            outf << " Апроксимация линейной функции  = "<< R;
            cout << " Апроксимация линейной функции = "<< R;
        }


    }

    void LangrangeInterpol() {
        vector<double> x;
        vector<double> y;

        ifstream inf("info_in.txt");
        ifstream infn("info_in_y.txt");
        ifstream langrange("langrange_int.txt");
        ofstream outf("info_out.txt");
        if (!inf) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else if (!infn) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else if (!langrange) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else {
            double xk;
            int n;
            inf >> n;
            langrange >> xk;
            double yk = 0;
            cout << endl;
            for (int i = 0; i < n; i++) {
                double k;
                double t;
                inf >> k;
                infn >> t;
                x.push_back(k);
                y.push_back(t);
            }
            for(int i=0;i<n;i++){
                double P = 1;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        P = P * ((xk - x[j]) / (x[i] - x[j]));
                    }
                }
                yk = yk + P * y[i];
            }
            outf <<" Полином Лагранжа = "<< yk;
            cout << " Полином Лагранжа = " << yk;
       }
    }

    void KWInterpol() {
        vector<double> x;
        vector<double> y;

        ifstream inf("info_in.txt");
        ifstream infn("info_in_y.txt");
        ifstream langrange("langrange_int.txt");
        ofstream outf("info_out.txt");
        if (!inf) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else if (!infn) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else if (!langrange) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else {
            double xk;
            int n;
            inf >> n;
            langrange >> xk;
            cout << endl;
            outf << "Кусково квадратическая интерполяция = ";
            double yk;
            for (int i = 0; i < n; i++) {
                double k;
                inf >> k;
                x.push_back(k);
            }
            for (int j = 0; j < n; j++) {
                double t;
                infn >> t;
                y.push_back(t);
            }
            
            for (int i = 0; i < n - 2; i++) {
                if (xk <= x[i + 1]) {
                    break;
                }
                
                yk = (xk - x[i + 1]) * (xk - x[i + 2]) * y[i] / ((x[i] - x[i + 1]) * (x[i] - x[i + 2])) +
                    (xk - x[i]) * (xk - x[i + 2]) * y[i + 1] / ((x[i + 1] - x[i]) * (x[i + 1] - x[i + 2])) +
                    (xk - x[i]) * (xk - x[i + 1]) * y[i + 2] / ((x[i + 2] - x[i]) * (x[i + 2] - x[i + 1]));

            }
            cout << " Кусочно-квдратическая интерполяция = " << yk;
            outf << " Кусочно-квдратическая интерполяция = " << yk;


        }
    }
    void LinInterol() {
        vector<double> x;
        vector<double> y;

        ifstream inf("info_in.txt");
        ifstream infn("info_in_y.txt");
        ifstream langrange("langrange_int.txt");
        ofstream outf("info_out.txt");
        if (!inf) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else if (!infn) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else if (!langrange) {
            cerr << "Не удалось считать файл! " << endl;
            exit(1);
        }
        else {
            double xk;
            int n;
            inf >> n;
            langrange >> xk;
            for (int i = 0; i < n; i++) {
                double k;
                inf >> k;
                x.push_back(k);
            }
            for (int j = 0; j < n; j++) {
                double t;
                infn >> t;
                y.push_back(t);
            }
            double yk;
            for (int i = 0; i < n - 1; i++) {
                if (xk <= x[i + 1]) {
                    break;
                }
                yk = (xk - x[i + 1]) * y[i] / (x[i] - x[i + 1]) +
                    (xk - x[i]) * y[i + 1] / (x[i + 1] - x[i]);

            }

            outf << " Кусочно-линейная интерполяция =  " << yk;
            cout << endl;
            cout << " Кусочно-линейная интерполяция =  " << yk << endl;
        }
    }


    };

    int main()
    {
        setlocale(LC_ALL, "Rus");
        cout << "\n \t\t\tМетоди апроксимации й интерполяции функций, ихняя оценка" << endl << "\n";
        cout << "\t\t\t 1. Вывести информацию про введенные данные." << endl;
        cout << "\t\t\t 2. Апроксимация линейной функции." << endl;
        cout << "\t\t\t 3. Интеполяционный полином Лагранжа." << endl;
        cout << "\t\t\t 4. Кусочно квадратическая интерполяция." << endl;
        cout << "\t\t\t 5. Кусочно  линейная интерполяция." << endl;
        cout << endl;
        cout << " Введите число нужного задания: ";
        int n;
        cin >> n;
        Approximation numb1;
        if (n == 1) {
            numb1.print_array();
        }
        else if (n == 2) {
            numb1.Aprox_linear_func();
        }
        else if (n == 3) {
            numb1.LangrangeInterpol();
        }
        else if (n == 4) {
            numb1.KWInterpol();
        }
        else if (n == 5) {
            numb1.LinInterol();
        }
        return 0;
    }
