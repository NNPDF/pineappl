/**
 * @file PineAPPL.hpp
 * @brief object oriented interface to PineAPPL
 */
#ifndef PineAPPL_HPP_
#define PineAPPL_HPP_

#include <LHAPDF/LHAPDF.h>
#include <cstdint>
#include <pineappl_capi.h>
#include <string>

/** @brief object oriented interface to PineAPPL  */
namespace PineAPPL {

/** @brief Key-value storage for passing optional information during grid
 * creation */
struct KeyVal {
  /** @brief underlying raw object */
  pineappl_keyval *raw;

  /** @brief constructor */
  KeyVal() : raw(pineappl_keyval_new()) {}

  KeyVal(const KeyVal &) = delete;

  KeyVal(KeyVal &&) = delete;

  KeyVal &operator=(const KeyVal &) = delete;

  KeyVal &operator=(KeyVal &&) = delete;

  /** @brief destructor */
  ~KeyVal() { pineappl_keyval_delete(this->raw); }

  /** @name setter */
  ///@{
  void set_double(const std::string &key, const double value) const {
    pineappl_keyval_set_double(this->raw, key.c_str(), value);
  }
  void set_bool(const std::string &key, const bool value) const {
    pineappl_keyval_set_bool(this->raw, key.c_str(), value);
  }
  void set_int(const std::string &key, const int value) const {
    pineappl_keyval_set_int(this->raw, key.c_str(), value);
  }
  void set_string(const std::string &key, const std::string &value) const {
    pineappl_keyval_set_string(this->raw, key.c_str(), value.c_str());
  }
  ///@}

  /** @name getter */
  ///@{
  double get_double(const std::string &key) const {
    return pineappl_keyval_double(this->raw, key.c_str());
  }
  bool get_bool(const std::string &key) const {
    return pineappl_keyval_bool(this->raw, key.c_str());
  }
  int get_int(const std::string &key) const {
    return pineappl_keyval_int(this->raw, key.c_str());
  }
  std::string get_string(const std::string &key) const {
    return pineappl_keyval_string(this->raw, key.c_str());
  }
  ///@}
};

/** @brief Entry in luminosity function  */
struct LumiEntry {
  /** @brief first parton id */
  std::int32_t pid1;

  /** @brief second parton id */
  std::int32_t pid2;

  /** @brief relative weight */
  double weight;
};

/** @brief Luminosity function */
struct Lumi {
  /** @brief underlying raw object */
  pineappl_lumi *raw;

  /** @brief constructor */
  Lumi() : raw(pineappl_lumi_new()) {}

  Lumi(const Lumi &) = delete;

  Lumi(Lumi &&) = delete;

  Lumi &operator=(const Lumi &) = delete;

  Lumi &operator=(Lumi &&) = delete;

  /** @brief destructor */
  ~Lumi() { pineappl_lumi_delete(this->raw); }

  /** @brief number of elements */
  std::size_t count() const { return pineappl_lumi_count(this->raw); }

  /**
   * @brief add a luminosity function
   * @param c luminosity function
   */
  void add(const std::vector<LumiEntry> &c) const {
    const std::size_t n = c.size();
    std::vector<std::int32_t> pids;
    std::vector<double> weights;
    for (const LumiEntry &l : c) {
      pids.push_back(l.pid1);
      pids.push_back(l.pid2);
      weights.push_back(l.weight);
    }
    pineappl_lumi_add(this->raw, n, pids.data(), weights.data());
  }

  /**
   * @brief Returns the number of combinations of the luminosity function `lumi`
   * for the specified entry.
   * @param entry position in lumi
   */
  std::size_t combinations(const std::size_t entry) const {
    return pineappl_lumi_combinations(this->raw, entry);
  }
};

/** @brief Coupling powers for each grid. */
struct Order {
  /** @brief Exponent of the strong coupling. */
  std::uint32_t alphas;

  /** @brief Exponent of the electromagnetic coupling. */
  std::uint32_t alpha;

  /** @brief Exponent of the logarithm of the scale factor of the renomalization
   * scale. */
  std::uint32_t logxir;

  /** @brief Exponent of the logarithm of the scale factor of the factorization
   * scale. */
  std::uint32_t logxif;
};

struct Grid {

  /** @brief underlying raw object */
  pineappl_grid *raw;

  /**
   * @brief constructor
   * @param lumi luminosity
   * @param orders orders
   * @param bin_limits bin limits
   * @param key_val additional informations
   */
  Grid(const Lumi &lumi, const std::vector<Order> &orders,
       const std::vector<double> &bin_limits, const KeyVal &key_val) {
    // cast orders
    const std::size_t n_orders = orders.size();
    std::vector<std::uint32_t> raw_orders;
    for (const Order &order : orders) {
      raw_orders.push_back(order.alphas);
      raw_orders.push_back(order.alpha);
      raw_orders.push_back(order.logxir);
      raw_orders.push_back(order.logxif);
    }
    this->raw = pineappl_grid_new(lumi.raw, n_orders, raw_orders.data(),
                                  bin_limits.size() - 1, bin_limits.data(),
                                  key_val.raw);
  }

  Grid() = delete;

  Grid(const Grid &) = delete;

  Grid(Grid &&) = delete;

  Grid &operator=(const Grid &) = delete;

  Grid &operator=(Grid &&) = delete;

  /** @brief destructor */
  ~Grid() { pineappl_grid_delete(this->raw); }

  /**
   * @brief Number of orders
   * @return number of orders
   */
  std::size_t order_count() const {
    return pineappl_grid_order_count(this->raw);
  }

  /**
   * @brief Number of bins
   * @return number of bins
   */
  std::size_t bin_count() const { return pineappl_grid_bin_count(this->raw); }

  /**
   * @brief Fill grid for the given parameters
   * @param x1 first momentum fraction
   * @param x2 second momentum fraction
   * @param q2 scale
   * @param order order index
   * @param observable observable value
   * @param lumi luminosity index
   * @param weight weight
   */
  void fill(const double x1, const double x2, const double q2,
            const std::size_t order, const double observable,
            const std::size_t lumi, const double weight) const {
    pineappl_grid_fill(this->raw, x1, x2, q2, order, observable, lumi, weight);
  }

  /**
   * @brief perform a convolution of the grid with PDFs
   * @param pdg_id hadron ID
   * @param pdf PDF
   * @param xi_ren renormalization scale variation
   * @param xi_fac factorization scale variation
   * @param order_mask order mask
   * @param lumi_mask luminosity mask
   * @return prediction for each bin
   */
  std::vector<double>
  convolute_with_one(const std::int32_t pdg_id, LHAPDF::PDF &pdf,
                     const double xi_ren = 1.0, const double xi_fac = 1.0,
                     const std::vector<bool> &order_mask = {},
                     const std::vector<bool> &lumi_mask = {}) const {
    // prepare LHAPDF stuff
    auto xfx = [](std::int32_t id, double x, double q2, void *pdf) {
      return static_cast<LHAPDF::PDF *>(pdf)->xfxQ2(id, x, q2);
    };
    auto alphas = [](double q2, void *pdf) {
      return static_cast<LHAPDF::PDF *>(pdf)->alphasQ2(q2);
    };
    // cast order_mask
    std::unique_ptr<bool[]> raw_order_mask;
    if (!order_mask.empty()) {
      raw_order_mask = std::unique_ptr<bool[]>(new bool[order_mask.size()]);
      std::copy(order_mask.begin(), order_mask.end(), &raw_order_mask[0]);
    }
    // cast lumi mask
    std::unique_ptr<bool[]> raw_lumi_mask;
    if (!lumi_mask.empty()) {
      raw_lumi_mask = std::unique_ptr<bool[]>(new bool[lumi_mask.size()]);
      std::copy(lumi_mask.begin(), lumi_mask.end(), &raw_lumi_mask[0]);
    }
    // do it!
    std::vector<double> results(this->bin_count());
    pineappl_grid_convolute_with_one(this->raw, pdg_id, xfx, alphas, &pdf,
                                     raw_order_mask.get(), raw_lumi_mask.get(),
                                     xi_ren, xi_fac, results.data());
    return results;
  }

  /**
   * @brief Write grid to file
   * @param filename file name
   */
  void write(const std::string &filename) const {
    pineappl_grid_write(this->raw, filename.c_str());
  }

  /**
   * @brief Set a metadata entry
   * @param key key
   * @param value value
   */
  void set_key_value(const std::string &key, const std::string &value) const {
    pineappl_grid_set_key_value(this->raw, key.c_str(), value.c_str());
  }

  /**
   * @brief Get a metadata entry
   * @param key key
   * @return value
   */
  std::string get_key_value(const std::string &key) const {
    auto *value = pineappl_grid_key_value(this->raw, key.c_str());
    std::string res(value);
    // delete the allocated object
    pineappl_string_delete(value);
    return res;
  }
};

} // namespace PineAPPL

#endif // PineAPPL_HPP_
