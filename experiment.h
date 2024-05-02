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
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
namespace macrodr {

using var::Var;
using var::Vector_Space;

class Time : public Var<Time, double> {};
class number_of_samples : public var::Constant<number_of_samples, double> {};

class ATP_concentration : public Var<ATP_concentration, double> {};

using ATP_step = var::Vector_Space<number_of_samples, ATP_concentration>;

using ATP_evoltype = std::variant<ATP_step, std::vector<ATP_step>>;

class ATP_evolution : public Var<ATP_evolution, std::vector<ATP_step>> {
public:
    auto size() const { return (*this)().size(); }
    
    bool is_zero() const {
        for (std::size_t i = 0; i < size(); ++i)
            if (get<ATP_concentration>((*this)[i])() != 0)
                return false;
        return true;
    }
    
    ATP_step &operator[](std::size_t i) { return (*this)()[i]; }
    
    ATP_step const &operator[](std::size_t i) const { return (*this)()[i]; }
    
    using Var<ATP_evolution, std::vector<ATP_step>>::Var;
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
    return put(f, sep, i_frac, i_step, time, v());
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
  std::string line;
  std::getline(f, line);
  std::stringstream ss(line);
  
  ss >> septr("i_step") >> sep >> septr("patch_current");
  if (!ss)
    return error_message("not even titles");
  std::getline(f, line);
  ss = std::stringstream(line);
  
  while (extract_double(ss >> i_step >> sep, current, separator[0])) {
      if (i_step != out().size())
          return error_message("i_step mismatch " + std::to_string(i_step) +
                               " size: " + std::to_string(out().size()));
      else
          out().push_back(Patch_current(current));
      std::getline(f, line);
      ss = std::stringstream(line);
  }
  return out;
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

        if ((i_frac > 0) && (i_frac != v.size() + 1))
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
      step().push_back(ATP_step(number_of_samples(std::size_t(v_ns)),
                                ATP_concentration(v_ATP)));
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
                            Recording_conditions const &r) {
    std::ofstream f(filename);
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    f << "i_step" << sep << "time" << sep << "i_sub_step" << sep
      << "number_of_samples" << sep << "ATP_concentration"
      << "\n";
    
    
    for (auto i = 0ul; i < r().size(); ++i) {
        Experiment_step const &s = r()[i];
        auto &t = get<Time>(s);
        auto const &a = get<ATP_evolution>(s);
        for (auto j = 0ul; j < a().size(); ++j)
            f << i << sep << t << sep << j << sep << get<number_of_samples>(a()[j])
              << sep << get<ATP_concentration>(a()[j]) << "\n";
    }
}


inline void save_experiment(const std::string filename, std::string sep,
                            Experiment const &e) {
    save_experiment(filename,sep,get<Recording_conditions>(e));
}

inline void save_fractioned_experiment(const std::string filename,
                                       std::string sep,
                                       std::vector<Experiment> const &e) {
  std::ofstream f(filename);
  f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  f << "i_frac" << sep << "i_step" << sep << "time" << sep << "i_sub_step"
    << sep << "number_of_samples" << sep << "ATP_concentration"
    << "\n";

  for (std::size_t k = 0; k < e.size(); ++k) {
    auto &r = get<Recording_conditions>(e[k]);
    auto fs = get<Frequency_of_Sampling>(e[k]);

    for (auto i = 0ul; i < r().size(); ++i) {
      Experiment_step const &s = r()[i];
      auto &t = get<Time>(s);
      auto const &av = get<ATP_evolution>(s)();
                     for (auto j = 0ul; j < av.size(); ++j)

                  f << k << sep << i << sep
                    << (j > 0 ? t() + get<number_of_samples>(av[j - 1])() / fs()
                              : t())
                    << sep << j << sep << get<number_of_samples>(av[j]) << sep
                    << get<ATP_concentration>(av[j]) << "\n";
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
            if (i_frac >= i_frac_prev) {
                if (as.size() == 0)
                    return error_message("mismatch here");
                else {
                    rec().emplace_back(Time(t_prev), as);
                    as.clear();
                }
                i_step_prev = i_step;
                i_sub_step_prev = std::numeric_limits<std::size_t>::max();
            }
        }
        if (i_frac_prev != i_frac) {
            
            if ((i_frac > 0) && (i_frac != e.size() + 1))
                return error_message("i_frac missmatch expected " +
                                     std::to_string(e.size() + 1) +
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
        } else
            return error_message("i_sub_step missmatch expected" +
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
    rec().emplace_back(Time(t_prev), as);
    e.emplace_back(rec, Frequency_of_Sampling(frequency_of_sampling),
                   initial_ATP_concentration(ATP_concentration(initial_ATP)));
    
    return true;
}

inline Maybe_error<bool>
load_experiment(const std::string filename, std::string separator,
                           double frequency_of_sampling, double initial_ATP,
                           
                           Experiment &e) {
    
    std::ifstream f(filename);
    if (!f)
        return error_message("cannot open file " + filename);
    std::string line;
    std::getline(f, line);
    std::stringstream ss(line);
    
    if (!(ss >> septr("i_step") >>
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
    std::size_t i_frac=0;
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
    
    while ((ss >> i_step >> septr(separator) >> t >>
            septr(separator) >> i_sub_step >> septr(separator) >>
            v_number_of_samples >> septr(separator) >> val)) {
        if (i_step_prev != i_step) {
            if ((i_step > 0) && (i_step != rec().size() + 1))
                return error_message("i_step missmatch expected " +
                                     std::to_string(rec().size() + 1) +
                                     " found:" + std::to_string(i_step));
            if (i_frac >= i_frac_prev) {
                if (as.size() == 0)
                    return error_message("mismatch here");
                else {
                    rec().emplace_back(Time(t_prev), as);
                    as.clear();
                }
                i_step_prev = i_step;
                i_sub_step_prev = std::numeric_limits<std::size_t>::max();
            }
        }
        if (i_frac_prev != i_frac) {
            
            if ((i_frac > 0) && (i_frac !=  1))
                return error_message("i_frac missmatch expected " +
                                     std::to_string( 1) +
                                     " found:" + std::to_string(i_frac));
            
            i_frac_prev = i_frac;
            rec().clear();
            i_step_prev = std::numeric_limits<std::size_t>::max();
        }
        
        if (i_sub_step_prev != i_sub_step) {
            if (i_sub_step != as.size())
                return error_message("i_sub_step missmatch expected" +
                                     std::to_string(as.size()) +
                                     " found:" + std::to_string(i_sub_step));
        } else
            return error_message("i_sub_step missmatch expected" +
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
    rec().emplace_back(Time(t_prev), as);
    e=Experiment(rec, Frequency_of_Sampling(frequency_of_sampling),
                   initial_ATP_concentration(ATP_concentration(initial_ATP)));
    
    return true;
}


// inline
//     std::pair<std::size_t, std::size_t>
// next_ATP_Pulse(const Recording_conditions& x, std::size_t pos)
// {
//     while (get<ATP_evolution>(x()[pos])
// }

// inline Recording_conditions
// Idealize_ATP_pulse(std::string save_name,Recording_conditions experiment) {
//     auto filename = save_name ;

//     // save_fractioned_experiment(filename + "_experiment.csv", ",", xs);
//     // save_fractioned_Recording(filename + "_recording.csv", ",", ys);
//     // return std::tuple(filename + "_experiment.csv", filename +
//     "_recording.csv",
//     //                   get<Frequency_of_Sampling>(experiment)(),
//     //                   get<initial_ATP_concentration>(experiment)()());
// }




struct Idealize_Pulses{
/**
 * nput vector of ATP_steps
| detect pulses (>threshold, in this case 0, but for experiment might be different)
| measure pulses area
| measure pulses maximum
-> area/maximum ->width
-> pulse center= (t50 decay +t50 rise)/2
-> pulse start= pulse center- half width
->pulse end = pulse center - half width

output: in the sample where the jump starts/ends, we divide the ATP_step in two:
0 until the t_stat
maxATP from there until t_stop
 * 
 * */




struct pulse_pos{
    std::size_t i_start;
    std::size_t i_end;
    
};

static auto measure_area_pulse(const Recording_conditions& x, pulse_pos p)
{
    
    double sum=0;
    auto i_run=p.i_start;
    while (i_run<p.i_end)
    {
        auto A=get<ATP_concentration>(get<ATP_evolution>(x()[i_run])()[0]);
        auto n=get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0]);
        sum+=A()*n();
        ++i_run;    
    }
    return sum;    
}

static auto measure_max_pulse(const Recording_conditions& x, pulse_pos p)
{
    double max=0;
    auto i_run=p.i_start;
    while (i_run<p.i_end)
    {
        auto A=get<ATP_concentration>(get<ATP_evolution>(x()[i_run])()[0]);
        if (A()>max)
            max=A();
        ++i_run;
    }
    return max;    
}


static auto find_t_50(const Recording_conditions& x, pulse_pos p, double max)
{
    bool found=false;
    double cum_samples=0;
    double t_50_up;
    double t_50_down;
    auto i_run=p.i_start;
    double A0=0;
    double A1;
    while (!found && i_run<p.i_end)
    {
        A1=get<ATP_concentration>(get<ATP_evolution>(x()[i_run])()[0])();
        if (A1>max/2.0){
            found=true;
            double dA=A1-A0;
            double nf=0.5-(A1-max/2.0)/dA;
            t_50_up=cum_samples+nf*get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0])();
            //cum_samples+=get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0])();
        }
        else{
        cum_samples+=get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0])();
        ++i_run;
        A0=A1;}
    }
    found=false;
    while (!found && i_run<p.i_end)
    {
        A1=get<ATP_concentration>(get<ATP_evolution>(x()[i_run])()[0])();
        if (A1==max){
            found=true;
          }
        else{
            cum_samples+=get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0])();
            ++i_run;
            }
    }
    found=false;
    while (!found && i_run<p.i_end)
    {
        A1=get<ATP_concentration>(get<ATP_evolution>(x()[i_run])()[0])();
        if (A1<max/2.0){
            found=true;
            double dA=A1-A0;
            double nf=0.5-(A1-max/2.0)/dA;
            t_50_down=cum_samples+nf*get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0])();
            //cum_samples+=get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0])();
        }
        else{
            cum_samples+=get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0])();
            ++i_run;
            A0=A1;
            t_50_down=cum_samples;
        }
    }
    
    
    return std::pair(t_50_up,t_50_down);
}


static bool starts_pulse(const ATP_evolution& x)
{
    return get<ATP_concentration>(x()[0])() >0;
}

static bool ends_pulse(const ATP_evolution& x)
{
    return get<ATP_concentration>(x()[0])() ==0;
}


static auto detect_Pulses(const Recording_conditions& x)
{
    assert([](auto x){
        for (auto& e: x())
            if (get<ATP_evolution>(e).size()>1)
                return false;
        return true;
    }(x));
    
    std::vector<pulse_pos> out;
    auto i_run=0ul;
    pulse_pos p;  
    while (i_run<x().size())
    {
        while (i_run<x().size()&&!starts_pulse(x()[i_run]))
        {
            ++i_run;
        }
        p.i_start=i_run;
        while (i_run<x().size()&&!ends_pulse(x()[i_run]))
        {
            ++i_run;
        }
        p.i_end=i_run;
        out.push_back(p);
    }
    return out;
}

static Recording_conditions& idealize_Pulse( Recording_conditions& x, pulse_pos pos)
{
    auto area=measure_area_pulse(x,pos);
    auto max=measure_max_pulse(x,pos);
    auto w= area/max;
    auto[t50u,t50d]=find_t_50(x,pos,max);
    auto tcenter=(t50u+t50d)/2;
    auto t_start=tcenter-w/2;
    auto t_end=tcenter+w/2;
    bool found=false;
    double cum_samples=0;
    auto i_run=pos.i_start;
    while (!found && i_run<pos.i_end)
    {
        auto n=get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0])();
        if (cum_samples+n>t_start){
            found=true;
            double n_0=(t_start-cum_samples);
            double n_1=n-n_0;
            auto a0=ATP_step(number_of_samples(n_0),ATP_concentration(0.0));
            auto a1=ATP_step(number_of_samples(n_1),ATP_concentration(max));
            get<ATP_evolution>(x()[i_run])()[0]=a0;
            assert(get<ATP_evolution>(x()[i_run])().size()==1);
            get<ATP_evolution>(x()[i_run])().push_back(a1);
            cum_samples+=n;
            ++i_run;
        }
        else{
            get<ATP_concentration>(get<ATP_evolution>(x()[i_run])()[0])()=0.0;
            cum_samples+=get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0])();
            ++i_run;
        }
    }
    found=false;
    while (!found && i_run<pos.i_end)
    {
        auto n=get<number_of_samples>(get<ATP_evolution>(x()[i_run])()[0])();
        if (cum_samples+n>t_end){
            found=true;
            double n_0=(t_end-cum_samples);
            double n_1=n-n_0;
            auto a0=ATP_step(number_of_samples(n_0),ATP_concentration(max));
            auto a1=ATP_step(number_of_samples(n_1),ATP_concentration(0.0));
            get<ATP_evolution>(x()[i_run])()[0]=a0;
            assert(get<ATP_evolution>(x()[i_run])().size()==1);
            get<ATP_evolution>(x()[i_run])().push_back(a1);
            cum_samples+=n;
            ++i_run;
        }
        else{
            get<ATP_concentration>(get<ATP_evolution>(x()[i_run])()[0])()=max;
            cum_samples+=n;
            ++i_run;
        }
    }
    while (i_run<pos.i_end)
    {
        get<ATP_concentration>(get<ATP_evolution>(x()[i_run])()[0])()=0.0;
        ++i_run;
    }
    return x;
  }
  
  
  
  auto operator()(const Recording_conditions& x)
{
    auto pulses=detect_Pulses(x);
    auto out=x;
    for (auto& e: pulses)
        out=idealize_Pulse(out,e);
    return out;
          
}
inline auto idealize_Experiment_Pulse(const Experiment& x) {
    auto out=x;
    
    get<Recording_conditions>(out)=Idealize_Pulses{}(get<Recording_conditions>(x));
    
    return out;
}

};




namespace cmd {
inline auto get_Experiment(
    std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt",
    double frequency_of_sampling = 50e3, double initial_ATP = 0) {
  using namespace macrodr;
    Experiment e;
  auto Maybe_exp=load_experiment(filename,",",frequency_of_sampling,initial_ATP,e);
    if (Maybe_exp)
  {
      return e;
  }
    else
  {
  auto [recording_conditions, recording] = macrodr::load_recording(filename);

  return Experiment(std::move(recording_conditions),
                    Frequency_of_Sampling(frequency_of_sampling),
                    initial_ATP_concentration(ATP_concentration(initial_ATP)));
    }

    }
    
    
    inline auto get_Experiment_file(
        std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt",
        double frequency_of_sampling = 50e3, double initial_ATP = 0) {
        return std::tuple(filename,frequency_of_sampling,initial_ATP);
    }
    
    
    inline void idealize_Experiment(
    std::string experiment ,
    std::string sep,
    std::string output) {
    using namespace macrodr;
    
    auto [recording_conditions, recording] = macrodr::load_recording(experiment);
    
    recording_conditions=Idealize_Pulses{}(recording_conditions);
    
    save_experiment(output,sep, recording_conditions);
 }


inline auto get_Observations(
    std::string filename = "../macro_dr/Moffatt_Hume_2007_ATP_time_7.txt") {
  auto [recording_conditions, recording] = macrodr::load_recording(filename);
  std::string out = filename.substr(0, filename.size() - 4) + "_recording.txt";
  save_Recording(out, ",", recording);
  return out;
}
using recording_type =
    typename return_type<std::decay_t<decltype(&get_Observations)>>::type;

using experiment_type =
    typename return_type<std::decay_t<decltype(&get_Experiment)>>::type;

using experiment_file_type =
    typename return_type<std::decay_t<decltype(&get_Experiment_file)>>::type;

} // namespace cmd

} // namespace macrodr

#endif // EXPERIMENT_H
