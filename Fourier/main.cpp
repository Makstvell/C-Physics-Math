#include <bits/stdc++.h>
using namespace std;

using cd = complex<double>;
const double PI = acos(-1.0);

// In-place iterative FFT (Cooley–Tukey, radix-2)
// invert=false -> FFT
// invert=true  -> IFFT (після чого ділимо на N)
void fft(vector<cd>& a, bool invert) {
    int n = (int)a.size();

    // 1) Bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }

    // 2) Layers: len = 2,4,8,...,n
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? +1 : -1);
        cd wlen(cos(ang), sin(ang)); // e^{+- i 2π/len}

        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i + j];
                cd v = a[i + j + len/2] * w;

                a[i + j] = u + v;
                a[i + j + len/2] = u - v;

                w *= wlen;
            }
        }
    }

    // 3) Normalize for inverse
    if (invert) {
        for (cd & x : a) x /= n;
    }
}

// Зручна обгортка: перетворити vector<double> у спектр
vector<cd> fft_real(const vector<double>& x) {
    vector<cd> a(x.begin(), x.end());
    // N має бути степінь 2 (або доповнюємо нулями)
    fft(a, false);
    return a;
}

int main() {
    // Приклад: 8 точок (N=2^k)
    vector<double> x = {1, 0, 0, 0, 0, 0, 0, 0};

    vector<cd> a(x.begin(), x.end());

    cout << "Init:\n";
    for (int k = 0; k < (int)a.size(); k++) {
        cout << k << ": " << a[k] << "\n";
    }


    fft(a, false);

    cout << "FFT:\n";
    for (int k = 0; k < (int)a.size(); k++) {
        cout << k << ": " << a[k] << "\n";
    }

    // IFFT повертає назад
    fft(a, true);
    cout << "\nBack (IFFT):\n";
    for (int n = 0; n < (int)a.size(); n++) {
        cout << n << ": " << a[n] << "\n";
    }

    return 0;
}
