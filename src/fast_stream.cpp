#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // Para std::nth_element

using namespace Rcpp;

// Estruturas internas

struct LiftingStep {
  std::string type;
  std::vector<double> coeffs;
  int start_idx;
};

// Classe que gerencia memoria e calculos com Zero-Allocation
class WaveletEngine {
public:
  // Buffers principais
  std::vector<double> buffer; // Janela deslizante (Ring Buffer)
  std::vector<double> dense_sig; // Sinal linearizado para processamento

  // Espacos de trabalho para LWT (separados para evitar alocacao)
  std::vector<double> work_even;
  std::vector<double> work_odd;

  // Espaco de rascunho para estatistica (MAD)
  // calcular mediana exige alterar a ordem dos dados.
  // Usamos este vetor "sujo" para nao estragar os dados reais.
  std::vector<double> scratch_mad;

  // Configuracao
  std::vector<LiftingStep> steps;
  double norm_approx, norm_detail;
  int levels;
  int window_size;
  int ext_mode;

  // Estado
  int head;
  int count;

  // Cache de Thresholds
  std::vector<double> current_lambdas;

  // Construtor
  WaveletEngine(
    List r_steps, NumericVector norm,
    int lvl, int w_size,
    int mode
    ) {
    window_size = w_size;
    levels = lvl;
    ext_mode = mode;

    norm_approx = norm[0];
    norm_detail = norm[1];

    // Parse dos passos
    int n_steps = r_steps.size();
    for(int i=0; i<n_steps; i++) {
      List s = r_steps[i];
      LiftingStep step;
      step.type = as<std::string>(s["type"]);
      step.coeffs = as<std::vector<double>>(s["coeffs"]);
      step.start_idx = s["start_idx"];
      steps.push_back(step);
    }

    // Alocacao única (Zero-Alloc durante execucao)
    buffer.resize(window_size, 0.0);
    dense_sig.resize(window_size);

    int half_size = (window_size + 1) / 2;
    work_even.resize(half_size);
    work_odd.resize(half_size);
    scratch_mad.resize(half_size); // Tamanho máximo necessario (d1)

    // Inicializa lambdas com 0
    current_lambdas.resize(levels, 0.0);

    head = 0;
    count = 0;
  }

  // Helper de acesso (padding virtual)
  inline double get_val(const std::vector<double>& x, int i, int n) {
    if (i >= 0 && i < n) return x[i];

    if (ext_mode == 3) return 0.0; // Zero

    if (ext_mode == 2) { // Periodic
      int idx = i % n;
      if (idx < 0) idx += n;
      return x[idx];
    }

    // Symmetric
    while (i < 0 || i >= n) {
      if (i < 0) i = -1 - i;
      else i = 2 * n - 1 - i;
    }
    // Clamp
    if (i < 0) i = 0;
    if (i >= n) i = n - 1;
    return x[i];
  }

  // Funcao interna: calcular estatistica e thresholds
  void update_thresholds(double alpha, double beta) {
    int n_d1 = work_odd.size();
    // work_odd contem d1 apos o forward transform

    // Calcular MAD (Median Absolute Deviation)
    // Copia para scratch para poder ordenar sem destruir d1
    for(int i=0; i<n_d1; i++) {
      scratch_mad[i] = std::abs(work_odd[i]);
    }

    // Algoritmo O(N): std::nth_element
    // Coloca a mediana na posicao correta sem ordenar tudo
    int mid = n_d1 / 2;
    std::nth_element(
      scratch_mad.begin(), scratch_mad.begin() + mid,
      scratch_mad.end()
      );
    double mad_val = scratch_mad[mid];

    // Estimar sigma
    double sigma = mad_val / 0.6745;

    if (sigma < 1e-15) {
      std::fill(current_lambdas.begin(), current_lambdas.end(), 0.0);
      return;
    }

    //  Calcular lambda nivel 1 (Eq. 9 Liu et al)
    // lambda = beta * sigma * sqrt(2 * log(N))
    double lambda_1 = beta * sigma *
      std::sqrt(2.0 * std::log((double)n_d1));
    current_lambdas[0] = lambda_1;

    //  Recursao para niveis superiores (se houver)
    for (int k = 1; k < levels; k++) {
      int lvl = k + 1;
      double prev = current_lambdas[k-1];
      // L_i = L_{i-1} * (i-1) / (i + alpha - 1)
      double factor = (double)(lvl - 1) / (double)(lvl + alpha - 1);
      current_lambdas[k] = prev * factor;
    }
  }

  // Processamento principal
  double push_and_process(
      double new_val, double alpha,
      double beta, std::string method,
      int update_freq, int step_iter
  ) {
    // Ingestao circular
    buffer[head] = new_val;
    head = (head + 1) % window_size;
    if (count < window_size) count++;

    // Bypass (warming phase)
    if (count < window_size) return new_val;

    // Linearizar (ring -> dense)
    for(int i=0; i<window_size; i++) {
      dense_sig[i] = buffer[(head + i) % window_size];
    }

    // Forward LWT (nivel 1 fixo neste loop otimizado)
    // Split
    int n = window_size;
    int n_even = (n + 1) / 2;
    int n_odd = n / 2;

    for(int i=0; i<n_even; i++) work_even[i] = dense_sig[2*i];
    for(int i=0; i<n_odd; i++)  work_odd[i]  = dense_sig[2*i+1];

    // Apply Steps (P/U)
    for(const auto& step : steps) {
      int k_filt = step.coeffs.size();
      if (step.type == "predict") {
        for(int i=0; i<n_odd; i++) {
          double sum = 0.0;
          for(int k=0; k<k_filt; k++) {
            sum += get_val(work_even, i + step.start_idx + k, n_even) *
              step.coeffs[k];
          }
          work_odd[i] -= sum;
        }
      } else {
        for(int i=0; i<n_even; i++) {
          double sum = 0.0;
          for(int k=0; k<k_filt; k++) {
            sum += get_val(work_odd, i + step.start_idx + k, n_odd) *
              step.coeffs[k];
          }
          work_even[i] += sum;
        }
      }
    }

    // Normaliza
    for(int i=0; i<n_even; i++) work_even[i] *= norm_approx;
    for(int i=0; i<n_odd; i++)  work_odd[i]  *= norm_detail;

    // Verifica se precisa atualizar estatisticas (Lazy Check)
    if (step_iter % update_freq == 0) {
      update_thresholds(alpha, beta);
    }

    // Aplica threshold (hard/soft/semisoft) usando o valor atual
    // (em cache ou novo)
    double lambda = current_lambdas[0]; // pegando nivel 1
    double lambda_sq = lambda * lambda;

    for(int i=0; i<n_odd; i++) {
      double v = work_odd[i];
      double av = std::abs(v);

      if (av < lambda) {
        work_odd[i] = 0.0;
      } else {
        if (method == "soft") {
          work_odd[i] = (v > 0) ? (av - lambda) : -(av - lambda);
        } else if (method == "semisoft") {
          double s = std::sqrt(v*v - lambda_sq);
          work_odd[i] = (v > 0) ? s : -s;
        }
        // hard: mantem valor
      }
    }

    // Inverse ILWT
    // Denormaliza
    for(int i=0; i<n_even; i++) work_even[i] /= norm_approx;
    for(int i=0; i<n_odd; i++)  work_odd[i]  /= norm_detail;

    // Reverse Steps
    for (int j = steps.size() - 1; j >= 0; j--) {
      const auto& step = steps[j];
      int k_filt = step.coeffs.size();

      if (step.type == "predict") {
        for(int i=0; i<n_odd; i++) {
          double sum = 0.0;
          for(int k=0; k<k_filt; k++) {
            sum += get_val(work_even, i + step.start_idx + k, n_even) *
              step.coeffs[k];
          }
          work_odd[i] += sum;
        }
      } else {
        for(int i=0; i<n_even; i++) {
          double sum = 0.0;
          for(int k=0; k<k_filt; k++) {
            sum += get_val(work_odd, i + step.start_idx + k, n_odd) *
              step.coeffs[k];
          }
          work_even[i] -= sum;
        }
      }
    }

    // Merge e retorna ultimo ponto (causalidade)
    // Indice original n-1.
    // Se n-1 par -> even[(n-1)/2].
    // Se n-1 impar -> odd[(n-1-1)/2] = odd[(n-2)/2].
    int last_idx = window_size - 1;
    if (last_idx % 2 == 0) {
      return work_even[last_idx / 2];
    } else {
      return work_odd[(last_idx - 1) / 2];
    }
  }
};

// EXPORTS PARA O R

//' Turbo causal batch
//' @keywords internal
// [[Rcpp::export]]
 NumericVector run_turbo_batch(
     NumericVector signal, List steps, NumericVector norm,
     int window_size, double alpha, double beta,
     std::string method, int ext_mode,
     int update_freq
 ) {
   int n = signal.size();
   NumericVector output(n);

   // Instancia Engine
   // Nota: fixamos levels=1 para esta implementacao turbo focada em
   // Denoising padrao.
   // Se quiser multi-level no futuro, precisa expandir o loop interno.
   WaveletEngine engine(steps, norm, 1, window_size, ext_mode);

   for(int i=0; i<n; i++) {
     // Passamos 'i' como step_iter para controlar o Lazy Update
     output[i] = engine.push_and_process(
       signal[i], alpha,
       beta, method, update_freq, i
     );
   }

   return output;
 }
