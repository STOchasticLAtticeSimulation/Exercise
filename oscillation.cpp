// 単振動を解いて t, x, xdot を "oscillation.dat" に出力する

using namespace std;

void RK4forHO(double &t, double &x, double &xdot, double dt); // t, x, xdot を渡すと単振動 EoM に従い dt だけ t, x, xdot を更新する

// ------------ パラメータ ----------------- //
const string filename = "oscillation.dat"; // 出力ファイル名
const double tf = 10; // 終了時刻
const double dt = 0.01; // 時間刻み
// ----------------------------------------- //

int main()
{}
