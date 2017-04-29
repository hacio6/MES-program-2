#include <iostream>
#include <fstream>

using namespace std;

void wypisz(double **macierz, double *wyniki, int wymiar);
bool gauss(int n, double **AB, double *X);

struct element {
	int Np;
	int WB;
	double K;
	double C;
	double Ro;
	double Rmax;
	double Rp;
	double Tot;
	double Alfa;
	double dR;
	double dTau;
	double TpTau;
	double *E, *W, *N1, *N2;
	double **Ke;
	double *Fe;
};
int main() {

	int nh = 9;
	int ne = nh - 1;
	int Np = 2;
	double Rmin, Rmax, AlfaAir, Tp, Tot, C, Ro, K, TauMax;
	double dR, a;
	double E[2], W[2], N1[2], N2[2], r[2];
	double x;
	int pom;
	
	/*-------------------------------------------------------------------------------------
	---------------------------WCZYTYWANIE DANYCH Z PLIKU--------------------------------
	-------------------------------------------------------------------------------------*/

	fstream p, wynik;
	p.open("dane.txt", ios::in | ios::out);
	if (p.good() == true) {
		p >> Rmin;
		p >> Rmax;
		p >> AlfaAir;
		p >> Tp;
		p >> Tot;
		p >> C;
		p >> Ro;
		p >> K;
		p >> TauMax;

		cout << "Dane:" << endl;
		cout << "Rmin - promien minimalny [m] : " << Rmin << endl;
		cout << "Rmax - promien maksymalny [m] : " << Rmax << endl;
		cout << "AlfaAir - wspolczynnik konwekcyjnej wymiany ciepla [W/(m^2*C)] : " << AlfaAir << endl;
		cout << "Tp - temperatura poczatkowa [K] : " << Tp << endl;
		cout << "Tot - temperatura otoczenia [K] : " << Tot << endl;
		cout << "C - cieplo wlasciwe [J/(kg*K)] : " << C << endl;
		cout << "Ro - gestosc [kg/m^3] : " << Ro << endl;
		cout << "K - wspolczynnik przewodzenia ciepla [W/(m*K)] : " << Rmin << endl;
		cout << "TauMax - czas procesu [s] : " << TauMax << endl;
	}
	double **gK = new double *[nh];
	for (int i = 0; i < nh; i++) {
		gK[i] = new double[nh];
	}
	double *gF = new double[nh];

	element *el = new element[nh];
	a = K / (C * Ro);
	W[0] = 1;
	W[1] = 1;
	E[0] = -0.5773502692;
	E[1] = 0.5773502692;
	N1[0] = 0.5 * (1 - E[0]);
	N1[1] = 0.5 * (1 - E[1]);
	N2[0] = 0.5 * (1 + E[0]);
	N2[1] = 0.5 * (1 + E[1]);
	dR = (Rmax - Rmin) / ne;
	

	/*-------------------------------------------------------------------------------------
	--------------------------------------MACIERZE---------------------------------------
	-------------------------------------------------------------------------------------*/
	for (int i = 0; i < nh; i++) {
		el[i].WB = 0;
	}
	el[ne-1].WB = 1;
	for (int i = 0; i < ne; i++)
	{
		el[i].Ke = new double *[ne];
		for (int j = 0; j < 2; j++)
		{
			el[i].Ke[j] = new double[2];
		}

	}
	for (int i = 0; i < ne; i++)
	{
		el[i].E = new double[2];
		el[i].W = new double[2];
		el[i].N1 = new double[2];
		el[i].N2 = new double[2];
		el[i].Fe = new double[2];
	}

	
	for (int dTau = 50; dTau <= TauMax; dTau += 50) {
		cout << endl;
		cout << "dTau [s] : " << dTau << endl;
		for (int i = 0; i < nh; i++) {
			for (int j = 0; j < nh; j++) {
				gK[i][j] = 0;
				gF[j] = 0;
			}
		}
		x = 0; 
		pom = 0;
		for (int i = 0; i < ne; i++)
		{
			
			r[0] = x;
			r[1] = x + dR;
			x = x + dR;
			el[i].Np = Np;
			el[i].K = K;
			el[i].C = C;
			el[i].Ro = Ro;
			el[i].Rmax = Rmax;
			el[i].Tot = Tot;
			el[i].N1[0] = N1[0];
			el[i].N1[1] = N1[1];
			el[i].N2[0] = N2[0];
			el[i].N2[1] = N2[1];
			el[i].dR = dR;
			el[i].W[0] = W[0];
			el[i].W[1] = W[1];
			el[i].dTau = dTau;
			el[i].Alfa = AlfaAir;
			el[i].Ke[0][0] = 0;
			el[i].Ke[0][1] = 0;
			el[i].Ke[1][0] = 0;
			el[i].Ke[1][1] = 0;
			el[i].Fe[0] = 0;
			el[i].Fe[1] = 0;
			el[i].Rp = 0;
			for (int j = 0; j < 2; j++) {
				el[i].Rp = el[i].N1[j] * r[0] + el[i].N2[j] * r[1];
				el[i].TpTau = el[i].N1[j] * Tp + el[i].N2[j] * Tp;
				el[i].Ke[0][0] = el[i].Ke[0][0] + el[i].K * el[i].Rp * el[i].W[j] / el[i].dR + el[i].C * el[i].Ro * el[i].dR * el[i].Rp * el[i].W[j] * el[i].N1[j] * el[i].N1[j] / el[i].dTau;
				el[i].Ke[0][1] = el[i].Ke[0][1] - el[i].K * el[i].Rp * el[i].W[j] / el[i].dR + el[i].C * el[i].Ro * el[i].dR * el[i].Rp * el[i].W[j] * el[i].N1[j] * el[i].N2[j] / el[i].dTau;
				el[i].Ke[1][0] = el[i].Ke[0][1];
				el[i].Ke[1][1] = el[i].Ke[1][1] + (el[i].K * el[i].Rp * el[i].W[j]) / el[i].dR + el[i].C * el[i].Ro * el[i].dR * el[i].Rp * el[i].W[j] * el[i].N2[j] * el[i].N2[j] / el[i].dTau;
				
				el[i].Fe[0] = el[i].Fe[0] - el[i].C * el[i].Ro * el[i].dR * el[i].TpTau * el[i].Rp * el[i].W[j] * el[i].N1[j] / el[i].dTau;
				el[i].Fe[1] = el[i].Fe[1] - el[i].C * el[i].Ro * el[i].dR * el[i].TpTau * el[i].Rp * el[i].W[j] * el[i].N2[j] / el[i].dTau;			
			}
			if (el[i].WB == 1) {
				el[i].Ke[1][1] = el[i].Ke[1][1] + 2 * el[i].Alfa*el[i].Rmax;
				el[i].Fe[1] = el[i].Fe[1] - (2 * el[i].Alfa * Rmax * Tot);
			}
			
			//---------------------- GENEROWANIE MACIERZY GLOBALNEJ --------------------------
			
			gF[i] = gF[i] + el[i].Fe[0];
			gF[i + 1] = gF[i + 1] + el[i].Fe[1];
			gK[pom][pom] = gK[pom][pom] + el[i].Ke[0][0];			
			gK[pom][pom + 1] = gK[pom][pom + 1] + el[i].Ke[0][1];			
			gK[pom + 1][pom] = gK[pom + 1][pom] + el[i].Ke[1][0];
			gK[pom + 1][pom + 1] = gK[pom + 1][pom + 1] + el[i].Ke[1][1];
			
			pom++;
		}
		wypisz(gK, gF, nh);

		double **HF;
		HF = new double *[nh];
		for (int i = 0; i < (nh); i++)
			HF[i] = new double[nh+1];

		for (int i = 0; i < (nh); i++)
			for (int j = 0; j <= (nh); j++)
				HF[i][j] = 0;


		for (int i = 0; i < (nh); i++)
		{
			for (int j = 0; j <= (nh); j++)
			{
				HF[i][j] = gK[i][j];
				HF[i][nh] = -gF[i];
			}
		}

		double *t;
		t = new double[nh];
		for (int i = 0; i < nh; i++) {
			t[i] = 0;
		}
		cout << "WYNIKI" << endl;
		wynik.open("wyniki.txt", ios::in | ios::out | ios::app);
		if (wynik.good() == true) {
			if (gauss(nh, HF, t)) {
				for (int i = 0; i < (nh); i++) {
					cout << "t" << i + 1 << " = " << t[i] << endl;
					wynik << t[i] << " ";
				} wynik << endl;
			}
			else
				cout << "COS NIE TAK!\n";
			wynik.close();
		}
		else cout << "Brak dostepu do pliku!" << endl;
		p.close();
		//gauss2(gK, gF, nh);
	}
	system("PAUSE");
	return 0;
}
void wypisz(double **macierz, double *wyniki, int wymiar) {
	cout << endl;
	for (int i = 0; i < wymiar; i++){
		for (int j = 0; j < wymiar; j++){
			cout.width(8);
			cout << macierz[i][j] << "  ";
			
		}
		cout << "  " << wyniki[i] << endl;
	}
	
	cout << endl;
	
}
bool gauss(int n, double **AB, double *X) // metoda algorytmu eliminacji Gaussa
{
	const double eps = 1e-12; // sta³a przybli¿enia zera
	int i, j, k;
	double m, s;

	// eliminacja wspó³czynników

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}


