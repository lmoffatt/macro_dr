#pragma once
#ifndef QMODEL_H
#define QMODEL_H
#include <map>
#include <string>

#include "Maybe_error.h"

namespace macrodr {
using logic::error_message;
using logic::Maybe_error;

class State_Model;

Maybe_error<State_Model>
to_State_Model(std::size_t t_number_of_states,
               std::map<std::pair<std::size_t, std::size_t>, std::string>
                   &t_transition_rates,
               std::map<std::pair<std::size_t, std::size_t>, std::string>
                   &t_agonist_transition_rates,
               std::map<std::size_t, std::string> t_conductances) ;


class State_Model {
  std::size_t m_number_of_states;

  std::map<std::pair<std::size_t, std::size_t>, std::string>
      m_transition_rates;
  std::map<std::pair<std::size_t, std::size_t>, std::string>
      m_agonist_transition_rates;
  std::map<std::size_t, std::string> m_conductances;

  State_Model(std::size_t t_number_of_states,
              std::map<std::pair<std::size_t, std::size_t>, std::string>
                  &t_transition_rates,
              std::map<std::pair<std::size_t, std::size_t>, std::string>
                  &t_agonist_transition_rates,
              std::map<std::size_t, std::string> t_conductances)
      : m_number_of_states{t_number_of_states},
        m_transition_rates{t_transition_rates},
        m_agonist_transition_rates{t_agonist_transition_rates},
        m_conductances{t_conductances} {}

public:
  friend Maybe_error<State_Model>
  to_State_Model(std::size_t t_number_of_states,
                 std::map<std::pair<std::size_t, std::size_t>, std::string>
                     &t_transition_rates,
                 std::map<std::pair<std::size_t, std::size_t>, std::string>
                     &t_agonist_transition_rates,
                 std::map<std::size_t, std::string> t_conductances) {
    for (auto &elem : t_transition_rates) {
      if (elem.first.first >= t_number_of_states)
        return error_message(
            "transition start state greater than number of "
            "states; number_of_states = " +
            std::to_string(t_number_of_states) +
            " start state= " + std::to_string(elem.first.first));
      else if (elem.first.second >= t_number_of_states)
        return error_message(
            "transition end state greater than number of "
            "states; number_of_states = " +
            std::to_string(t_number_of_states) +
            " end state= " + std::to_string(elem.first.second));
      else if (elem.first.second == t_number_of_states)
        return error_message(
            "transition start state same as end state;  start state: " +
            std::to_string(elem.first.first) +
            " end state= " + std::to_string(elem.first.second));
    }
    //
    for (auto &elem : t_agonist_transition_rates) {
      if (elem.first.first >= t_number_of_states)
        return error_message(
            "agonist transition start state greater than number of "
            "states; number_of_states = " +
            std::to_string(t_number_of_states) +
            " start state= " + std::to_string(elem.first.first));
      else if (elem.first.second >= t_number_of_states)
        return error_message(
            "agonist transition end state greater than number of "
            "states; number_of_states = " +
            std::to_string(t_number_of_states) +
            " end state= " + std::to_string(elem.first.second));
      else if (elem.first.second == t_number_of_states)
        return error_message(
            "agonist transition start state same as end state;  start state: " +
            std::to_string(elem.first.first) +
            " end state= " + std::to_string(elem.first.second));
    }

    std::size_t i_cond = 0;
    for (auto &elem : t_conductances) {
      if (elem.first >= t_number_of_states)
        return error_message(
            "state conductance number greater than number of states:"
            " proposed= " + std::to_string(elem.first) +
            " number_of_states = " + std::to_string(t_number_of_states));

      if (elem.first != i_cond)
        return error_message("state conductance skipped: current is" +
                             std::to_string(i_cond) +
                             " proposed= " + std::to_string(elem.first));
      ++i_cond;
    }
    if (i_cond != t_number_of_states)
      return error_message("state conductance missing: number of proposed states=" +
                            std::to_string(i_cond) +
                            " number_of_states = " + std::to_string(t_number_of_states));

    return State_Model( t_number_of_states,
                       t_transition_rates,
                       t_agonist_transition_rates,
                       t_conductances);


}
};

} // namespace macrodr

#endif // QMODEL_H
