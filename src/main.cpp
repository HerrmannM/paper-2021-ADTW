/** ADTW Demonstration application only
 *  ADTW Authors 2021
 */

#include "utils.hpp"
#include "tseries/tseries.hpp"
#include "tseries/readers/tsreader/tsreader.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <variant>
#include <mutex>
#include <limits>
#include <thread>
#include <queue>
#include <vector>
#include <chrono>

using namespace std;
namespace fs = std::filesystem;

constexpr auto INF = numeric_limits<double>::infinity();

[[nodiscard]] inline double gamma(double a, double b){
  double d = a-b;
  return d*d;
}

[[nodiscard]] inline double adtw(
  const double* lines, size_t nblines,
  const double* cols, size_t nbcols,
  const double weight,
  const double cutoff
) {

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // Create a new tighter upper bounds (most commonly uadtw in the code).
  // First, take the "next float" after "cutoff" to deal with numerical instability.
  // Then, subtract the cost of the last alignment.
  const double ub = nextafter(cutoff, INF)-gamma(lines[nblines-1], cols[nbcols-1]);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // Double buffer allocation, no initialisation required (border condition manage in the code).
  // Base indices for the 'c'urrent row and the 'p'revious row.
  auto buffers = std::unique_ptr<double[]>(new double[nbcols*2]);
  size_t c{0}, p{nbcols};

  // Line & column counters
  size_t i{0}, j{0};

  // Cost accumulator. Also uadtw as the "left neighbour".
  double cost;

  // EAP variables: track where to start the next line, and the position of the previous pruning point.
  // Must be init to 0: index 0 is the next starting index and also the "previous pruning point"
  size_t next_start{0}, prev_pp{0};

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // Initialisation of the first line.
  {
    const double l0 = lines[0];
    // Fist cell is a special case.
    // Check against the original upper bound dealing with the case where we have both series of length 1.
    cost = gamma(l0, cols[0]);
    if (cost>cutoff) { return INF; }
    buffers[c+0] = cost;
    // All other cells. Checking against "ub" is OK as the only case where the last cell of this line is the
    // last alignment is taken are just above (1==nblines==nbcols, and we have nblines >= nbcols).
    size_t curr_pp = 1;
    for (j = 1; j==curr_pp && j<nbcols; ++j) {
      cost = cost+gamma(l0, cols[j])+weight; // Left: penalty
      buffers[c+j] = cost;
      if (cost<=ub) { ++curr_pp; }
    }
    ++i;
    prev_pp = curr_pp;
  }

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // Main loop
  for (; i<nblines; ++i) {
    // --- --- --- Swap and variables init
    std::swap(c, p);
    const double li = lines[i];
    size_t curr_pp = next_start; // Next pruning point init at the start of the line
    j = next_start;
    // --- --- --- Stage 0: Special case for the first column. Can only look up (border on the left)
    {
      cost = buffers[p+j]+gamma(li, cols[j])+weight; // Top: penalty
      buffers[c+j] = cost;
      if (cost<=ub) { curr_pp = j+1; } else { ++next_start; }
      ++j;
    }
    // --- --- --- Stage 1: Up to the previous pruning point while advancing next_start: diag and top
    for (; j==next_start && j<prev_pp; ++j) {
      const auto d = gamma(li, cols[j]);
      cost = std::min(
        d+buffers[p+j-1],         // Diag: no penalty
        d+buffers[p+j]+weight     // Top: penalty
      );
      buffers[c+j] = cost;
      if (cost<=ub) { curr_pp = j+1; } else { ++next_start; }
    }
    // --- --- --- Stage 2: Up to the previous pruning point without advancing next_start: left, diag and top
    for (; j<prev_pp; ++j) {
      const auto d = gamma(li, cols[j]);
      cost = min(d+cost+weight,   // Left: penalty
        d+buffers[p+j-1],         // Diag: no penalty
        d+buffers[p+j]+weight);   // Top: penalty
      buffers[c+j] = cost;
      if (cost<=ub) { curr_pp = j+1; }
    }
    // --- --- --- Stage 3: At the previous pruning point. Check if we are within bounds.
    if (j<nbcols) { // If so, two cases.
      const auto d = gamma(li, cols[j]);
      if (j==next_start) { // Case 1: Advancing next start: only diag (no penalty)
        cost = buffers[p+j-1]+d;
        buffers[c+j] = cost;
        if (cost<=ub) { curr_pp = j+1; }
        else {
          // Special case if we are on the last alignment: return the actual cost if we are <= cutoff
          if (i==nblines-1 && j==nbcols-1 && cost<=cutoff) { return cost; }
          else { return INF; }
        }
      } else { // Case 2: Not advancing next start: possible path in previous cells: left (penalty) and diag.
        cost = std::min(d+cost+weight, d+buffers[p+j-1]);
        buffers[c+j] = cost;
        if (cost<=ub) { curr_pp = j+1; }
      }
      ++j;
    } else { // Previous pruning point is out of bound: exit if we extended next start up to here.
      if (j==next_start) {
        // But only if we are above the original UB
        // Else set the next starting point to the last valid column
        if (cost>cutoff) { return INF; }
        else { next_start = nbcols-1; }
      }
    }
    // --- --- --- Stage 4: After the previous pruning point: only prev.
    // Go on while we advance the curr_pp; if it did not advance, the rest of the line is guaranteed to be > ub.
    for (; j==curr_pp && j<nbcols; ++j) {
      const auto d = gamma(li, cols[j]);
      cost = cost+d+weight; // Left: penalty
      buffers[c+j] = cost;
      if (cost<=ub) { ++curr_pp; }
    }
    // --- --- ---
    prev_pp = curr_pp;
  } // End of main loop for(;i<nblines;++i)

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // Finalisation
  // Check for last alignment (i==nblines implied, Stage 4 implies j<=nbcols). Cost must be <= original bound.
  if (j==nbcols && cost<=cutoff) { return cost; } else { return INF; }
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
// Launch task in parallel
// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
class ParTasks {
public:

  using task_t = std::function<void()>;
  using taskgen_t = std::function<std::optional<task_t>()>;

  ParTasks() = default;

  /// Non thread safe! Add all the task before calling "execute"
  void push_task(task_t func) {
    tasklist.push(std::move(func));
  }

  /// Template version
  template<class F, class... Args>
  void push_task(F&& f, Args&& ... args) {
    tasklist.emplace(std::move(std::bind(f, args...)));
  }

  /// Blocking call
  void execute(int nbthreads) {
    if (nbthreads<=1) {
      while (!tasklist.empty()) {
        auto task = std::move(tasklist.front());
        tasklist.pop();
        task();
      }
    } else {
      threads.reserve(nbthreads);
      for (int i = 0; i<nbthreads; ++i) { threads.emplace_back([this]() { run_thread(); }); }
      // Wait for all threads to stop
      for (auto& thread : threads) { thread.join(); }
      threads.clear();
    }
  }

  /// Blocking call
  void execute(size_t nbthreads, size_t nbtask) {
    if (nbthreads<=1) {
      while (!tasklist.empty()) {
        auto task = std::move(tasklist.front());
        tasklist.pop();
        task();
      }
    } else {
      threads.reserve(nbthreads);
      for (size_t i = 0; i<nbthreads; ++i) { threads.emplace_back([this, nbtask]() { run_thread(nbtask); }); }
      // Wait for all threads to stop
      for (auto& thread : threads) { thread.join(); }
      threads.clear();
    }
  }

  /// Blocking call using a task generator
  void execute(size_t nbthread, taskgen_t tgenerator) {
    // --- --- --- 1 thread
    if (nbthread<=1) {
      auto ntask = tgenerator();
      while (ntask.has_value()) {
        auto task = ntask.value();
        task();
        ntask = tgenerator();
      }
    }
      // --- --- --- Multi thread
    else {
      threads.reserve(nbthread);
      for (size_t i = 0; i<nbthread; ++i) {
        threads.emplace_back([this, &tgenerator]() { run_thread_generator(tgenerator); });
      }
      // Wait for all threads to stop
      for (auto& thread : threads) { thread.join(); }
      threads.clear();
    }
  }

private:

  std::mutex mtx;
  std::vector<std::thread> threads;
  std::queue<task_t> tasklist;

  void run_thread() {
    mtx.lock();
    while (!tasklist.empty()) {
      auto task = std::move(tasklist.front());
      tasklist.pop();
      mtx.unlock();
      task();
      mtx.lock();
    }
    mtx.unlock();
  }

  void run_thread(size_t nbtask) {
    if (nbtask<=1) { run_thread(); }
    else {
      std::vector<task_t> tasks;
      tasks.reserve(nbtask);
      mtx.lock();
      while (!tasklist.empty()) {
        while (!tasklist.empty() && tasks.size()<nbtask) {
          tasks.emplace_back(std::move(tasklist.front()));
          tasklist.pop();
        }
        mtx.unlock();
        for (auto& t:tasks) { t(); }
        tasks.clear();
        mtx.lock();
      }
      mtx.unlock();
    }
  }

  void run_thread_generator(taskgen_t& tgenerator) {
    mtx.lock();
    auto ntask = tgenerator();
    mtx.unlock();
    while (ntask.has_value()) {
      auto task = ntask.value();
      task();
      {
        std::lock_guard lg(mtx);
        ntask = tgenerator();
      }
    }
  }

};






// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
// Main application
// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

void print_usage(ostream& out){
  out << "ADTW NN1 classifier." << endl;
  out << "Usage:" << endl;
  out << "  ./path/exec <ucr folder> <dataset name> <penalty> <threads>" << endl;
}

[[noreturn]]void print_error_exit(const string& msg, int exit_code) {
  std::cerr << "Error: " << msg << endl;
  print_usage(std::cerr);
  exit(exit_code);
}

[[noreturn]] void arg_error() {
  print_error_exit("argument parsing", 1);
}

// --- --- --- Argument checking/parsing
struct Args {
  fs::path ucr_path;
  string dsname;
  double penalty{};
  size_t nbthreads{};
};

Args read_args(int argc, char** argv) {
  vector<string> argList(argv, argv+argc);
  Args res;
  int i = 1;
  // ---
  if (i<argc) { res.ucr_path = fs::path(argList[i]); ++i; } else { arg_error(); }
  // ---
  if (i<argc) { res.dsname = argList[i]; ++i; } else { arg_error(); }
  // ---
  if (i<argc) { res.penalty = stod(argList[i]); ++i; } else { arg_error(); }
  // ---
  if (i<argc) { res.nbthreads = stoi(argList[i]); ++i; } else { arg_error(); }
  // ---
  return res;
}

struct ProgressMonitor {

  size_t total;

  explicit ProgressMonitor(size_t max){
    total = max;
  }

  // --- --- --- Print progress
  void print_progress(std::ostream& out, size_t nbdone){
    if(nbdone>0) {
      const size_t vprev = (nbdone-1)*100/total;
      const size_t vnow = nbdone*100/total;
      const size_t vnow_tenth = vnow/10;
      const size_t vprev_tenth = vprev/10;
      if (vprev<vnow) {
        if (vprev_tenth<vnow_tenth) { out << vnow_tenth*10 << "% "; } else { out << "."; }
        std::flush(out);
      }
    }
  }
};



std::ostream& operator<<(std::ostream& os, std::chrono::nanoseconds ns)
{
  using namespace std::chrono;
  using days = duration<int, std::ratio<86400>>;
  auto d = duration_cast<days>(ns);
  ns -= d;
  auto h = duration_cast<hours>(ns);
  ns -= h;
  auto m = duration_cast<minutes>(ns);
  ns -= m;
  auto s = duration_cast<seconds>(ns);
  ns -= s;

  std::optional<int> fs_count;
  switch (os.precision()) {
    case 9: fs_count = ns.count();
      break;
    case 6: fs_count = duration_cast<microseconds>(ns).count();
      break;
    case 3: fs_count = duration_cast<milliseconds>(ns).count();
      break;
  }

  char fill = os.fill('0');
  if (d.count())
    os << d.count() << "d ";
  if (d.count() || h.count())
    os << std::setw(2) << h.count() << ":";
  if (d.count() || h.count() || m.count())
    os << std::setw(d.count() || h.count() ? 2 : 1) << m.count() << ":";
  os << std::setw(d.count() || h.count() || m.count() ? 2 : 1) << s.count();
  if (fs_count.has_value())
    os << "." << std::setw(os.precision()) << fs_count.value();
  if (!d.count() && !h.count() && !m.count())
    os << "s";

  os.fill(fill);
  return os;
}



int main(int argc, char** argv) {

  // --- Args
  Args args = read_args(argc, argv);
  fs::path dir_path = args.ucr_path/args.dsname;
  fs::path train_path = dir_path/(args.dsname+"_TRAIN.ts");
  fs::path test_path = dir_path/(args.dsname+"_TEST.ts");

  // --- Train
  // --- --- Loading
  cout << "Loading Train: " << train_path << "..." << flush;
  auto itrain = ifstream(train_path);
  auto rtrain = TSReader::read(itrain);
  cout << " done." << endl;
  if (rtrain.index()==0) { print_error_exit(get<0>(rtrain), 2); }
  // --- --- Checking
  const auto& train = get<1>(rtrain);
  cout << "TRAIN: " << train.series.size() << " x ";
  if (train.has_equallength()) { cout << train.shortest_length; }
  else { cout << "[" << train.shortest_length << ".." << train.longest_length << "]"; }
  cout << endl;
  if(train.has_missings()){ print_error_exit("Missing data detected - aborting", 3); }

  // --- Test
  // --- --- Loading
  cout << "Loading Test: " << test_path << "..." << flush;
  auto itest = ifstream(test_path);
  auto rtest = TSReader::read(itest);
  cout << " done." << endl;
  if (rtest.index()==0) { print_error_exit(get<0>(rtest), 2); }
  // --- --- Checking
  const auto& test = get<1>(rtest);
  cout << "TRAIN: " << test.series.size() << " x ";
  if (test.has_equallength()) { cout << test.shortest_length; }
  else { cout << "[" << test.shortest_length << ".." << test.longest_length << "]"; }
  cout << endl;
  if(test.has_missings()){ print_error_exit("Missing data detected - aborting", 3); }

  // --- NN1
  std::mutex mutex;
  size_t nb_correct{0};
  size_t nb_done{0};
  ProgressMonitor pm(train.series.size()*test.series.size());

  // --- --- --- NN1 task
  auto nn1task = [penalty=args.penalty, &mutex, &pm, &nb_done, &train, &test](size_t idxTest) -> string {
    double bsf = INF;
    size_t bcandidate = train.series.size()+1;
    const auto& s_test = test.series[idxTest];
    // --- --- --- NN1 loop
    for (size_t idxTrain{0}; idxTrain<train.series.size(); idxTrain++) {
      const auto& s_train = train.series[idxTrain];
      double res = adtw(s_train.data(), s_train.length(), s_test.data(), s_test.length(), penalty, bsf);
      if (res<bsf) {
        bsf = res;
        bcandidate = idxTrain;
      }
      std::lock_guard lg(mutex);
      nb_done++;
      pm.print_progress(cout, nb_done);
    }
    // --- --- ---
    return train.series[bcandidate].label().value();
  };

  // --- --- --- Task generator
  auto nn1task_gen = [&mutex, &test, &nn1task, &nb_correct, idxTest = (size_t)0]() mutable {
    if (idxTest<test.series.size()) {
      ParTasks::task_t task = [&test, idxTest, &nn1task, &mutex, &nb_correct]() {
      string l = nn1task(idxTest);
      if (l==test.series[idxTest].label().value()) {
          std::lock_guard lg(mutex);
          nb_correct++;
        }
      };
      idxTest++;
      return std::optional<ParTasks::task_t>(task);
    } else { return std::optional<ParTasks::task_t>(); }
  };

  // --- --- ---
  ParTasks p;
  auto start = std::chrono::steady_clock::now();
  p.execute(args.nbthreads, nn1task_gen);
  auto stop = std::chrono::steady_clock::now();
  auto duration = stop-start;

  // --- Print results
  cout << endl;
  cout << "NN1 classification done in " << duration << endl;
  cout << "correct: " << nb_correct << "/" << test.series.size() << endl;
  cout << "accuracy: " << ((double)nb_correct)/((double)test.series.size()) << endl;
  return 0;

}
