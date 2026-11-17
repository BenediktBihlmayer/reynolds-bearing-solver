#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>

constexpr double PI = 3.14159265358979323846;

int main()
{
    // 1. Parameter in SI-Einheiten
    const double R = 0.026979;   // m
    const double B = 0.0182;     // m
    const double C = 0.000032;   // m
    const double eps = 0.66;
    const double omega = 200.0; // rad/s
    const double U = omega * R;
    const double eta = 0.010115; // Pa*s

    // Engstelle bei 315 Grad beziehungsweise -45 Grad
    const double phi_min_deg = 315.0;
    const double phi_min_rad = phi_min_deg * PI / 180.0;

    // Phasenverschiebung gegenueber dem Standardminimum bei PI
    const double delta_phi = phi_min_rad - PI;

    // Verschiebung fuer die Gitterverdichtung
    const double target_xi = phi_min_deg / 360.0;
    const double xi_shift = target_xi - 0.5;

    // 2. Gitter
    const int nx = 100;
    const int ny = 50;
    const double dz = B / (ny - 1);
    const double dz2 = dz * dz;

    std::vector<std::vector<double>> p(
        nx, std::vector<double>(ny, 0.0)
    );
    std::vector<double> phi(nx);
    std::vector<double> h(nx);
    std::vector<double> h_cubed(nx);

    const double alpha = 0.8;

    // Offset, damit das Gitter bei phi = 0 beginnt
    const double offset =
        alpha * std::sin(2.0 * PI * (0.0 - xi_shift));

    for (int i = 0; i < nx; ++i) {
        const double xi = static_cast<double>(i) / nx;

        phi[i] = 2.0 * PI * xi
               + alpha * std::sin(2.0 * PI * (xi - xi_shift))
               - offset;

        h[i] = C * (1.0 + eps * std::cos(phi[i] - delta_phi));
        h_cubed[i] = h[i] * h[i] * h[i];
    }

    // 3. Gauss-Seidel-Iteration
    const double tolerance = 1e-5;
    double max_error = 1.0;
    int iteration = 0;
    const int max_iterations = 100000;

    const double rhs_const = -6.0 * eta * U * C * eps / R;

    std::cout << "Starte Simulation ...\n";

    while (max_error > tolerance && iteration < max_iterations) {
        max_error = 0.0;

        for (int i = 0; i < nx; ++i) {
            // Periodische Randbedingungen in Umfangsrichtung
            const int i_prev = (i == 0) ? nx - 1 : i - 1;
            const int i_next = (i == nx - 1) ? 0 : i + 1;

            const double phi_next =
                (i == nx - 1) ? phi[0] + 2.0 * PI : phi[i_next];

            const double phi_prev =
                (i == 0) ? phi[nx - 1] - 2.0 * PI : phi[i_prev];

            const double dx_east = R * (phi_next - phi[i]);
            const double dx_west = R * (phi[i] - phi_prev);
            const double dx_avg = 0.5 * (dx_east + dx_west);

            const double h_east =
                0.5 * (h_cubed[i] + h_cubed[i_next]);
            const double h_west =
                0.5 * (h_cubed[i] + h_cubed[i_prev]);
            const double h_center = h_cubed[i];

            const double RHS =
                rhs_const * std::sin(phi[i] - delta_phi);

            const double A_E = h_east / (dx_east * dx_avg);
            const double A_W = h_west / (dx_west * dx_avg);
            const double A_N = h_center / dz2;
            const double A_S = h_center / dz2;

            for (int j = 0; j < ny; ++j) {
                // Axiale Neumann-Randbedingungen: dp/dz = 0
                const int j_prev = (j == 0) ? 1 : j - 1;
                const int j_next = (j == ny - 1) ? ny - 2 : j + 1;

                const double p_old = p[i][j];

                p[i][j] = (
                    A_E * p[i_next][j]
                    + A_W * p[i_prev][j]
                    + A_N * p[i][j_next]
                    + A_S * p[i][j_prev]
                    - RHS
                ) / (A_E + A_W + A_N + A_S);

                const double error = std::abs(p[i][j] - p_old);

                if (error > max_error) {
                    max_error = error;
                }
            }
        }

        ++iteration;

        if (iteration % 1000 == 0) {
            std::cout << "Iteration " << iteration
                      << ", max. Aenderung: " << max_error << "\n";
        }
    }

    if (max_error <= tolerance) {
        std::cout << "Konvergiert nach " << iteration
                  << " Iterationen.\n";
    } else {
        std::cerr << "Iterationslimit erreicht. Max. Aenderung: "
                  << max_error << "\n";
    }

    // 4. Lagerkraefte durch numerische Integration
    double F_X = 0.0;
    double F_Y = 0.0;

    for (int i = 0; i < nx; ++i) {
        const int i_prev = (i == 0) ? nx - 1 : i - 1;
        const int i_next = (i == nx - 1) ? 0 : i + 1;

        const double phi_next =
            (i == nx - 1) ? phi[0] + 2.0 * PI : phi[i_next];

        const double phi_prev =
            (i == 0) ? phi[nx - 1] - 2.0 * PI : phi[i_prev];

        const double dphi = 0.5 * (phi_next - phi_prev);

        for (int j = 0; j < ny; ++j) {
            const double weight_z =
                (j == 0 || j == ny - 1) ? 0.5 : 1.0;

            const double dA = R * dphi * dz * weight_z;

            F_X += -p[i][j] * std::cos(phi[i]) * dA;
            F_Y += -p[i][j] * std::sin(phi[i]) * dA;
        }
    }

    const double F_res = std::sqrt(F_X * F_X + F_Y * F_Y);
    const double gamma_rad = std::atan2(F_Y, F_X);
    const double gamma_deg = gamma_rad * 180.0 / PI;

    std::cout << "\nBerechnete Lagerkraefte:\n"
              << "F_X = " << F_X << " N\n"
              << "F_Y = " << F_Y << " N\n"
              << "F_res = " << F_res << " N\n"
              << "Schwenkwinkel gamma = " << gamma_deg
              << " Grad (" << gamma_rad << " rad)\n\n";

    // 5. CSV-Export der Druckverteilung
    std::cout << "Schreibe Druckdaten in CSV-Datei ...\n";
    std::ofstream outfile("gleitlager_druck.csv");

    if (!outfile.is_open()) {
        std::cerr << "Fehler: Druckdatei konnte nicht erstellt werden.\n";
        return 1;
    }

    for (int i = 0; i <= nx; ++i) {
        double phi_val;
        int i_idx;

        if (i == nx) {
            phi_val = 2.0 * PI;
            i_idx = 0;
        } else {
            phi_val = phi[i];
            i_idx = i;
        }

        for (int j = 0; j < ny; ++j) {
            const double z_val =
                static_cast<double>(j) / (ny - 1) * B;

            outfile << phi_val << ","
                    << z_val << ","
                    << p[i_idx][j] << "\n";
        }
    }
    outfile.close();

    // 6. CSV-Export der Kraefte
    std::cout << "Schreibe Kraefte in CSV-Datei ...\n";
    std::ofstream outfile_forces("gleitlager_kraefte.csv");

    if (!outfile_forces.is_open()) {
        std::cerr << "Fehler: Kraftdatei konnte nicht erstellt werden.\n";
        return 1;
    }

    outfile_forces
        << "F_X_N,F_Y_N,F_res_N,Winkel_gamma_rad,Winkel_gamma_deg\n"
        << F_X << "," << F_Y << "," << F_res << ","
        << gamma_rad << "," << gamma_deg << "\n";

    outfile_forces.close();

    // 7. Nachbearbeitung: Nulldurchgang auf PI verschieben
    std::cout << "\nBerechne Verschiebung fuer das Diagramm ...\n";

    const int j_mid = ny / 2;
    int max_idx = 0;
    double max_p = p[0][j_mid];

    for (int i = 0; i < nx; ++i) {
        if (p[i][j_mid] > max_p) {
            max_p = p[i][j_mid];
            max_idx = i;
        }
    }

    // Ab dem Maximum vorwaerts bis zum ersten Druck <= 0 suchen
    int null_idx = max_idx;

    for (int i = 0; i < nx; ++i) {
        const int current_idx = (max_idx + i) % nx;

        if (p[current_idx][j_mid] <= 0.0) {
            null_idx = current_idx;
            break;
        }
    }

    const double aktueller_null_winkel = phi[null_idx];
    const double shift = PI - aktueller_null_winkel;

    struct DataPoint {
        double phi_neu;
        double z;
        double p_val;
    };

    std::vector<DataPoint> plot_data;
    plot_data.reserve(nx * ny);

    for (int i = 0; i < nx; ++i) {
        const double phi_neu =
            std::fmod(phi[i] + shift + 4.0 * PI, 2.0 * PI);

        for (int j = 0; j < ny; ++j) {
            const double z_val =
                static_cast<double>(j) / (ny - 1) * B;

            plot_data.push_back({phi_neu, z_val, p[i][j]});
        }
    }

    // Strikte Sortierung nach Winkel und anschliessend z
    std::sort(
        plot_data.begin(),
        plot_data.end(),
        [](const DataPoint& a, const DataPoint& b) {
            if (a.phi_neu != b.phi_neu) {
                return a.phi_neu < b.phi_neu;
            }
            return a.z < b.z;
        }
    );

    std::ofstream outfile_plot("gleitlager_druck_verschoben.csv");

    if (!outfile_plot.is_open()) {
        std::cerr << "Fehler: Verschobene Druckdatei konnte "
                     "nicht erstellt werden.\n";
        return 1;
    }

    outfile_plot << "Winkel_neu_rad,Breite_z_m,Druck_p_Pa\n";

    for (const auto& point : plot_data) {
        outfile_plot << point.phi_neu << ","
                     << point.z << ","
                     << point.p_val << "\n";
    }

    // Originaler Abschluss des Diagramms.
    // ACHTUNG: Nach der Winkelverschiebung entspricht p[0][j]
    // im Allgemeinen nicht dem Druck am neuen Winkel 2*PI.
    for (int j = 0; j < ny; ++j) {
        const double z_val =
            static_cast<double>(j) / (ny - 1) * B;

        outfile_plot << 2.0 * PI << ","
                     << z_val << ","
                     << p[0][j] << "\n";
    }

    outfile_plot.close();

    std::cout << "Verschobene Druckverteilung gespeichert unter "
                 "'gleitlager_druck_verschoben.csv'.\n";

    return 0;
}
