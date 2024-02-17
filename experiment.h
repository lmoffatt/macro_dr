#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "general_output_operator.h"
#include "maybe_error.h"
#include "parallel_tempering.h"
// #include "qmodel.h"
#include "variables.h"
#include <cstddef>
#include <fstream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
namespace macrodr {

using var::Var;
using var::Vector_Space;

class Time : public Var<Time, double> {};
class number_of_samples : public var::Constant<number_of_samples, double> {};

class ATP_concentration : public Var<ATP_concentration, double> {};

using ATP_step = var::Vector_Space<number_of_samples, ATP_concentration>;

using ATP_evoltype = std::variant<ATP_step, std::vector<ATP_step>>;

class ATP_evolution
    : public Var<ATP_evolution, std::variant<ATP_step, std::vector<ATP_step>>> {
  using Var<ATP_evolution, std::variant<ATP_step, std::vector<ATP_step>>>::Var;
};

inline std::ostream &put(std::ostream &f, std::string sep, std::size_t i_frac,
                         std::size_t i_step, double time,
                         std::vector<ATP_step> const &v) {
  for (std::size_t i = 0; i < v.size(); ++i)
    f << i_frac << sep << i_step << sep << time << sep
      << i_step + (i + 0.5) / v.size() << sep << get<number_of_samples>(v[i])
      << sep << get<ATP_concentration>(v[i]) << "\n";
  return f;
}
inline std::ostream &put(std::ostream &f, std::string sep, std::size_t i_frac,
                         std::size_t i_step, double time, ATP_step const &x) {
  f << i_frac << sep << i_step << sep << time << sep << i_step + 0.5 << sep
    << get<number_of_samples>(x) << sep << get<ATP_concentration>(x) << "\n";
  return f;
}

inline std::ostream &put(std::ostream &f, std::string sep, std::size_t i_frac,
                         std::size_t i_step, double time,
                         ATP_evolution const &v) {
  return std::visit(
      [&f, sep, i_frac, i_step, time](auto &e) -> decltype(auto) {
        return put(f, sep, i_frac, i_step, time, e);
      },
      v());
}

class Patch_current : public Var<Patch_current, double> {};

using Experiment_step = var::Vector_Space<Time, ATP_evolution>;

class Recording : public Var<Recording, std::vector<Patch_current>> {};

class Recording_conditions
    : public Var<Recording_conditions, std::vector<Experiment_step>> {};

class Frequency_of_Sampling : public var::Var<Frequency_of_Sampling, double> {};

class initial_ATP_concentration
    : public Var<initial_ATP_concentration, ATP_concentration> {
    using Var<initial_ATP_concentration, ATP_concentration>::Var;
};

using Experiment =
    var::Vector_Space<Recording_conditions, Frequency_of_Sampling,
                      initial_ATP_concentration>;

} // namespace macrodr

template <class Parameter> class save_Parameter;
namespace macrodr {
inline auto &extract_double(std::istream &is, double &r, char sep) {
  std::stringstream ss;
    is.get(*ss.rdbuf(), sep);
  auto s = ss.str();
  if ((s == "nan") || (s == "NAN") || (s == "NaN"))
    r = std::numeric_limits<double>::quiet_NaN();
  else {
    try {
      r = std::stod(s);
    } catch (...) {
      is.setstate(std::ios::failbit);
    }
  }
  return is;
}

inline auto &extract_double(std::istream &is, double &r) {
    std::string s;
    is >> s;
    if ((s == "nan") || (s == "NAN") || (s == "NaN"))
        r = std::numeric_limits<double>::quiet_NaN();
    else {
        try {
            r = std::stod(s);
        } catch (...) {
            is.setstate(std::ios::failbit);
        }
    }
    return is;
}

template <class Parameter>
void report_model(save_Parameter<Parameter> &s,
                  std::vector<Experiment> const &e) {
  std::ofstream f(s.fname + "_experiment.csv");
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  f << "frequency_of_sampling\n" << get<Frequency_of_Sampling>(e[0])() << "\n";
  f << "initial_ATP_concentration\n"
    << get<initial_ATP_concentration>(e[0])() << "\n";

  f << "i_frac" << s.sep << "i_step" << s.sep << "time" << s.sep << "i_step_f"
    << s.sep << "number_of_samples" << s.sep << "ATP"
    << "\n";

  for (auto i_frac = 0ul; i_frac < e.size(); ++i_frac) {
    auto r = get<Recording_conditions>(e[i_frac]);
    for (auto i_step = 0ul; i_step < r().size(); ++i_step) {
      auto sa = r()[i_step];
      auto t = get<Time>(sa)();
      auto ev = get<ATP_evolution>(sa);
      put(f, s.sep, i_frac, i_step, t, ev);
    }
  }
}

template <class Parameter>
void report_model(save_Parameter<Parameter> &s,
                  std::vector<Recording> const &e) {
  std::ofstream f(s.fname + "_ifrac_recording.csv");
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  f << "i_frac" << s.sep << "i_step" << s.sep << "patch_current"
    << "\n";

  for (auto i_frac = 0ul; i_frac < e.size(); ++i_frac) {
    for (auto i_step = 0ul; i_step < e[i_frac]().size(); ++i_step) {
      f << i_frac << s.sep << i_step << s.sep << e[i_frac]()[i_step]() << "\n";
    }
  }
}

template <class Parameter>
Maybe_error<Recording> load_recording(save_Parameter<Parameter> &s) {
  std::ifstream f(s.fname + "_recording.csv");
  Recording out;
  std::size_t i_step;
  double current;
  f >> septr("i_step") >> s.sep >> septr("patch_current") >> septr("\n");
  if (!f)
    return error_message("not even titles");
  else {
    while (extract_double(f >> i_step >> s.sep, current)) {
      if (i_step != out().size())
        return error_message("i_step mismatch " + std::to_string(i_step) +
                             " size: " + std::to_string(out().size()));
      else
        out().push_back(Patch_current(current));
      std::string line;
      std::getline(f, line);
    }
    return out;
  }
}

inline Maybe_error<Recording> load_Recording(std::string filename,
                                             const std::string separator) {
  std::ifstream f(filename);
  Recording out;
  std::size_t i_step;
  auto sep = septr(separator);
  double current;
  f >> septr("i_step") >> sep >> septr("patch_current") >> septr("\n");
  if (!f)
    return error_message("not even titles");
  else {
    while (extract_double(f >> i_step >> sep, current, separator[0])) {
      if (i_step != out().size())
        return error_message("i_step mismatch " + std::to_string(i_step) +
                             " size: " + std::to_string(out().size()));
      else
        out().push_back(Patch_current(current));
      std::string line;
      std::getline(f, line);
    }
    return out;
  }
}

template <class Parameter>
void report_model(save_Parameter<Parameter> &s, Recording const &e) {
  std::ofstream f(s.fname + "_recording.csv");
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  f << "i_step" << s.sep << "patch_current"
    << "\n";

  for (auto i_step = 0ul; i_step < e().size(); ++i_step) {
    f << i_step << s.sep << e()[i_step]() << "\n";
  }
}

inline void save_Recording(std::string const &fname,
                           std::string const &separator, Recording const &e) {
  std::ofstream f(fname);
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  f << "i_step" << separator << "patch_current"
    << "\n";

  for (auto i_step = 0ul; i_step < e().size(); ++i_step) {
    f << i_step << separator << e()[i_step]() << "\n";
  }
}

inline void save_fractioned_Recording(std::string const &fname,
                                      std::string const &separator,
                                      std::vector<Recording> const &v) {
  std::ofstream f(fname);
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  f << "i_frac" << separator << "i_step" << separator << "patch_current"
    << "\n";
  for (auto i_frac = 0ul; i_frac < v.size(); ++i_frac) {
    {
      auto &e = v[i_frac];
      for (auto i_step = 0ul; i_step < e().size(); ++i_step) {
        f << i_frac << separator << i_step << separator << e()[i_step]()
          << "\n";
      }
    }
  }
}

inline Maybe_error<bool> load_fractioned_Recording(std::string const &fname,
                                                   std::string const &separator,
                                                   std::vector<Recording> &v) {

  std::ifstream f(fname);
  if (!f)
    return error_message("cannot open file " + fname);
  std::string line;
  std::getline(f, line);
  std::stringstream ss(line);

  if (!(ss >> septr("i_frac") >> septr(separator) >> septr("i_step") >>
        septr(separator) >> septr("patch_current")))
    return error_message("titles are wrong : expected  i_step:" + separator +
                         "patch_current; found:" + line);
  else {
    std::getline(f, line);
    ss = std::stringstream(line);
    std::size_t i_frac;
    std::size_t i_frac_prev = std::numeric_limits<std::size_t>::max();
    std::size_t i_step;
    std::size_t i_step_prev = std::numeric_limits<std::size_t>::max();
    Recording rec;
    double val;
    while (extract_double(ss >> i_frac >> septr(separator) >> i_step >>
                              septr(separator),
                          val, separator[0])) {
      if (i_frac_prev != i_frac) {
            
            if ((i_frac>0)&&(i_frac != v.size()+1))
          return error_message("i_frac missmatch expected" +
                               std::to_string(v.size()) +
                               " found:" + std::to_string(i_frac));
        else {
          if (i_frac > 0)
            v.push_back(rec);
          i_frac_prev = i_frac;
          rec().clear();
        }
      }

      if (i_step_prev != i_step) {
        if (i_step != rec().size())
          return error_message("i_step missmatch expected" +
                               std::to_string(rec().size()) +
                               " found:" + std::to_string(i_step));
        else {
          rec().push_back(Patch_current(val));
          i_step_prev = i_step;
        }
      }
      std::getline(f, line);
      ss = std::stringstream(line);
    }
    v.push_back(rec);
    
    return true;
  }
}

inline Maybe_error<bool> load_Recording_Data(std::string const &fname,
                                             std::string const &separator,
                                             Recording &e) {
  std::ifstream f(fname);
  if (!f)
    return error_message("cannot open file " + fname);
  else {
    std::string line;
    std::getline(f, line);
    std::stringstream ss(line);

    if (!(ss >> septr("i_step") >> septr(separator) >> septr("patch_current")))
      return error_message("titles are wrong : expected  i_step:" + separator +
                           "patch_current; found:" + line);
    else {
      std::getline(f, line);
      ss = std::stringstream(line);
      std::size_t i_step;
      std::size_t i_step_prev = std::numeric_limits<std::size_t>::max();

      double val;
      while (
          extract_double(ss >> i_step >> septr(separator), val, separator[0])) {
        if (i_step_prev != i_step) {
          if (i_step != e().size())
            return error_message("i_step missmatch expected" +
                                 std::to_string(e().size()) +
                                 " found:" + std::to_string(i_step));
          else {
            e().push_back(Patch_current(val));
            i_step_prev = i_step;
          }
        }
        std::getline(f, line);
        ss = std::stringstream(line);
      }
    }
    return true;
  }
}

inline std::tuple<Recording_conditions, Recording>
load_recording(const std::string filename) {
  std::ifstream f(filename);
  std::vector<Experiment_step> out0;
  std::vector<Patch_current> out1;

  std::string line;
  std::getline(f, line);
  while (std::getline(f, line)) {
    if (!line.empty()) {
      std::stringstream ss(line);
      double v_time, v_ns, v_ATP, v_current;
      extract_double(ss >> v_time >> v_ns >> v_ATP, v_current);
      ATP_evolution step;
      step() = ATP_step(number_of_samples(std::size_t(v_ns)),
                        ATP_concentration(v_ATP));
      out0.emplace_back(Time(v_time / 1000), std::move(step));
      out1.emplace_back(Patch_current(v_current));
    }
  }
  return std::tuple(Recording_conditions(out0), Recording(out1));
}

inline std::tuple<Experiment, Recording>
load_experiment(const std::string filename, double frequency_of_sampling,
                double initial_ATP) {
  auto [v_recording_conditions, v_recording] = load_recording(filename);

  return {Experiment(std::move(v_recording_conditions),
                     Frequency_of_Sampling(frequency_of_sampling),
                     initial_ATP_concentration(ATP_concentration(initial_ATP))),
          std::move(v_recording)};
}

inline void save_experiment(const std::string filename, std::string sep,
                            Experiment const &e) {
  std::ofstream f(filename);
  f << "i_step" << sep << "time" << sep << "i_sub_step" << sep
    << "number_of_samples" << sep << "ATP_concentration"
    << "\n";

  auto &r = get<Recording_conditions>(e);

  for (auto i = 0ul; i < r().size(); ++i) {
    Experiment_step const &s = r()[i];
    auto &t = get<Time>(s);
    auto const &a = get<ATP_evolution>(s);
    std::visit(overloaded(
                   [&t, &f, i, &sep](const ATP_step &a0) {
                     f << i << sep << t << sep << 0 << sep
                       << get<number_of_samples>(a0) << sep
                       << get<ATP_concentration>(a0) << "\n";
                   },
                   [&t, &f, i, &sep](const std::vector<ATP_step> &av) {
                     for (auto j = 0ul; j < av.size(); ++j)
                       f << i << sep << t << sep << j << sep
                         << get<number_of_samples>(av[j]) << sep
                         << get<ATP_concentration>(av[j]) << "\n";
                   }),
               a());
  }
}

inline void save_fractioned_experiment(const std::string filename,
                                       std::string sep,
                                       std::vector<Experiment> const &e) {
  std::ofstream f(filename);
  f << "i_frac" << sep << "i_step" << sep << "time" << sep << "i_sub_step"
    << sep << "number_of_samples" << sep << "ATP_concentration"
    << "\n";

  for (std::size_t k = 0; k < e.size(); ++k) {
    auto &r = get<Recording_conditions>(e[k]);
    auto fs = get<Frequency_of_Sampling>(e[k]);

    for (auto i = 0ul; i < r().size(); ++i) {
      Experiment_step const &s = r()[i];
      auto &t = get<Time>(s);
      auto const &a = get<ATP_evolution>(s);
      std::visit(
          overloaded(
              [&t, &f, k, i, &sep](const ATP_step &a0) {
                f << k << sep << i << sep << t() << sep << 0 << sep
                  << get<number_of_samples>(a0) << sep
                  << get<ATP_concentration>(a0) << "\n";
              },
              [&t, &f, i, k, &fs, &sep](const std::vector<ATP_step> &av) {
                for (auto j = 0ul; j < av.size(); ++j)

                  f << k << sep << i << sep
                    << (j > 0 ? t() + get<number_of_samples>(av[j - 1])() / fs()
                              : t())
                    << sep << j << sep << get<number_of_samples>(av[j]) << sep
                    << get<ATP_concentration>(av[j]) << "\n";
              }),
          a());
    }
  }
}

inline Maybe_error<bool>
load_fractioned_experiment(const std::string filename, std::string separator,
                           double frequency_of_sampling, double initial_ATP,

                           std::vector<Experiment> &e) {

  std::ifstream f(filename);
  if (!f)
    return error_message("cannot open file " + filename);
  std::string line;
  std::getline(f, line);
  std::stringstream ss(line);

  if (!(ss >> septr("i_frac") >> septr(separator) >> septr("i_step") >>
        septr(separator) >> septr("time") >> septr(separator) >>
        septr("i_sub_step") >> septr(separator) >> septr("number_of_samples") >>
        septr(separator) >> septr("ATP_concentration")))
    return error_message("titles are wrong : expected  "
                         "i_frac" +
                         separator + "i_step" + separator + "time" + separator +
                         "i_sub_step" + separator + "number_of_samples" +
                         separator +
                         "ATP_concentration"
                         "; found:" +
                         line);

  std::getline(f, line);
  ss = std::stringstream(line);
  std::size_t i_frac;
  std::size_t i_frac_prev = std::numeric_limits<std::size_t>::max();
  std::size_t i_step;
  std::size_t i_step_prev = std::numeric_limits<std::size_t>::max();
  double t;
  double t_prev = std::numeric_limits<double>::max();
  std::size_t i_sub_step;
  std::size_t i_sub_step_prev = std::numeric_limits<std::size_t>::max();
  double v_number_of_samples;
  Recording_conditions rec;
  std::vector<ATP_step> as;
  double val;

  while ((ss >> i_frac >> septr(separator) >> i_step >> septr(separator) >> t >>
          septr(separator) >> i_sub_step >> septr(separator) >>
          v_number_of_samples >> septr(separator) >> val)) {
      if (i_step_prev != i_step) {
          if ((i_step > 0) && (i_step != rec().size() + 1))
              return error_message("i_step missmatch expected " +
                                   std::to_string(rec().size() + 1) +
                                   " found:" + std::to_string(i_step));
         if  (i_frac>=i_frac_prev)
              {
                  if (as.size()==0)
                      return error_message("mismatch here");
                  else  {
                  if (as.size() == 1)
                      rec().emplace_back(Time(t_prev), as[0]);
                  else
                      rec().emplace_back(Time(t_prev), as);
                  as.clear();
              }
              i_step_prev = i_step;
              i_sub_step_prev = std::numeric_limits<std::size_t>::max();
           }
    }
    if (i_frac_prev != i_frac) {
          
          if ((i_frac>0)&&(i_frac != e.size()+1))
            return error_message("i_frac missmatch expected " +
                             std::to_string(e.size()+1) +
                             " found:" + std::to_string(i_frac));
      
        if (i_frac > 0)
              e.emplace_back(
                rec, Frequency_of_Sampling(frequency_of_sampling),
                initial_ATP_concentration(ATP_concentration(initial_ATP)));
        i_frac_prev = i_frac;
        rec().clear();
        i_step_prev = std::numeric_limits<std::size_t>::max();
      
    }

    if (i_sub_step_prev != i_sub_step) {
        if (i_sub_step != as.size())
        return error_message("i_sub_step missmatch expected" +
                             std::to_string(as.size()) +
                             " found:" + std::to_string(i_sub_step));
    }
    else return error_message("i_sub_step missmatch expected" +
                             std::to_string(as.size()) +
                             " found:" + std::to_string(i_sub_step));
    as.push_back(ATP_step(number_of_samples(v_number_of_samples),
                          ATP_concentration(val)));
    i_sub_step_prev = i_sub_step;
    if (i_sub_step_prev == 0)
      t_prev = t;
    std::getline(f, line);
    ss = std::stringstream(line);
  }
  return true;
}

namespace cmd {
inline auto get_Experiment(
    std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt",
    double frequency_of_sampling = 50e3, double initial_ATP = 0) {
    using namespace macrodr;
    auto [recording_conditions, recording] = macrodr::load_recording(filename);
    
    return Experiment(std::move(recording_conditions),
                      Frequency_of_Sampling(frequency_of_sampling),
                      initial_ATP_concentration(ATP_concentration(initial_ATP)));
}
inline auto get_Observations(
    std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt") {
    auto [recording_conditions, recording] = macrodr::load_recording(filename);
    std::string out = filename.substr(0, filename.size() - 4) + "_recording.txt";
    save_Recording(out, ",", recording);
    return out;
}

using experiment_type =
    typename return_type<std::decay_t<decltype(&get_Experiment)>>::type;

} // namespace cmd

} // namespace macrodr

#endif // EXPERIMENT_H
