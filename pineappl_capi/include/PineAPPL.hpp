/**
 * @brief Object oriented interface to PineAPPL.
 */
#ifndef PineAPPL_HPP_
#define PineAPPL_HPP_

#include <LHAPDF/LHAPDF.h>
#include <pineappl_capi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

/** @brief Object oriented interface to PineAPPL.  */
namespace PineAPPL {

/** @brief Key-value storage for passing optional information during grid
 * creation. */
struct KeyVal {
  /** @brief Underlying raw object. */
  pineappl_keyval *raw;

  /** @brief Constructor. */
  KeyVal() : raw(pineappl_keyval_new()) {}
  KeyVal(const KeyVal &) = delete;
  KeyVal(KeyVal &&) = delete;

  KeyVal &operator=(const KeyVal &) = delete;
  KeyVal &operator=(KeyVal &&) = delete;

  /** @brief Destructor. */
  ~KeyVal() { pineappl_keyval_delete(this->raw); }

  /** @name Setter. */
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

  /** @name Getter. */
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

/** @brief Entry in luminosity function.  */
struct LumiEntry {
  /** @brief First parton id. */
  std::int32_t pid1;

  /** @brief Second parton id. */
  std::int32_t pid2;

  /** @brief Relative weight. */
  double weight;
};

/** @brief Luminosity function. */
struct Lumi {
  /** @brief Underlying raw object. */
  pineappl_lumi *raw;

  /** @brief Constructor. */
  Lumi() : raw(pineappl_lumi_new()) {}
  Lumi(const Lumi &) = delete;
  Lumi(Lumi &&) = delete;

  Lumi &operator=(const Lumi &) = delete;
  Lumi &operator=(Lumi &&) = delete;

  /** @brief Destructor. */
  ~Lumi() { pineappl_lumi_delete(this->raw); }

  /** @brief Number of elements. */
  std::size_t count() const { return pineappl_lumi_count(this->raw); }

  /**
   * @brief Add a luminosity function.
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
   * @brief Returns the number of combinations of the luminosity function
   * `lumi` for the specified entry.
   * @param entry position in lumi
   */
  std::size_t combinations(const std::size_t entry) const {
    return pineappl_lumi_combinations(this->raw, entry);
  }
};

/** @brief Entry in sub-channel function.  */
struct SubChannelEntry {
  /** @brief First parton id. */
  std::vector<std::pair<std::vector<int>, double>> entry;
};

/** @brief Entry in channels function.  */
struct ChannelsEntry {
  /** @brief First parton id. */
  std::vector<SubChannelEntry> channels_entry;
};

/** @brief Channels function. */
struct Channels {
  /** @brief Underlying raw object. */
  pineappl_channels *raw;

  /** @brief Constructor. */
  Channels() : raw(pineappl_channels_new()) {}
  Channels(const Channels &) = delete;
  Channels(Channels &&) = delete;

  Channels &operator=(const Channels &) = delete;
  Channels &operator=(Channels &&) = delete;

  /** @brief Destructor. */
  ~Channels() { pineappl_channels_delete(this->raw); }

  /** @brief Number of elements. */
  std::size_t count() const { return pineappl_channels_count(this->raw); }

  /**
   * @brief Add a channel function.
   * @param c channel function
   */
  void add(const ChannelsEntry &c) const {
    const std::size_t combinations = c.channels_entry.size();
    const std::size_t nb_convolutions =
        c.channels_entry[0].entry[0].first.size();

    std::vector<std::int32_t> pids;
    std::vector<double> weights;
    for (const SubChannelEntry &s : c.channels_entry) {
      for (const std::pair<std::vector<int>, double> &m : s.entry) {
        for (const int &pid : m.first) {
          pids.push_back(pid);
        }
        weights.push_back(m.second);
      }
    }
    pineappl_channels_add(this->raw, combinations, nb_convolutions, pids.data(),
                          weights.data());
  }

  /**
   * @brief Returns the number of combinations of the channels function
   * `channels` for the specified entry.
   * @param entry position in channels
   */
  std::size_t combinations(const std::size_t entry) const {
    return pineappl_channels_combinations(this->raw, entry);
  }
};

/** @brief Coupling powers for each grid. */
struct Order {
  /** @brief Exponent of the strong coupling. */
  std::uint32_t alphas;

  /** @brief Exponent of the electromagnetic coupling. */
  std::uint32_t alpha;

  /** @brief Exponent of the logarithm of the scale factor of the
   * renomalization scale. */
  std::uint32_t logxir;

  /** @brief Exponent of the logarithm of the scale factor of the factorization
   * scale. */
  std::uint32_t logxif;
};

/** @brief Extension of `Order` to accommodate for v1 representation. */
struct OrderV1 {
  /** @brief Exponent of the strong coupling. */
  std::uint8_t alphas;

  /** @brief Exponent of the electromagnetic coupling. */
  std::uint8_t alpha;

  /** @brief Exponent of the logarithm of the scale factor of the
   * renomalization scale. */
  std::uint8_t logxir;

  /** @brief Exponent of the logarithm of the scale factor of the
   * factorization scale. */
  std::uint8_t logxif;

  /** @brief Exponent of the logarithm of the scale factor of the
   * fragmentation scale. */
  std::uint8_t logxia;
};

/** @brief Base Grid struct that contains the common data in v0 and v1. */
struct GridBase {
  /** @brief Underlying raw object. */
  pineappl_grid *raw;

  /** @brief Constructor (protected to avoid direct instantiation). */
 protected:
  GridBase(pineappl_grid *grid) : raw(grid) {}

  /** @brief Deleted copy/move semantics. */
  GridBase() = delete;
  GridBase(const GridBase &) = delete;
  GridBase(GridBase &&) = delete;

  GridBase &operator=(const GridBase &) = delete;
  GridBase &operator=(GridBase &&) = delete;

 public:
  /** @brief Destructor. */
  virtual ~GridBase() { pineappl_grid_delete(this->raw); }

  /**
   * @brief Number of orders.
   * @return number of orders
   */
  std::size_t order_count() const {
    return pineappl_grid_order_count(this->raw);
  }

  /**
   * @brief Number of bins.
   * @return number of bins
   */
  std::size_t bin_count() const { return pineappl_grid_bin_count(this->raw); }

  /**
   * @brief Write grid to file.
   * @param filename file name
   */
  void write(const std::string &filename) const {
    pineappl_grid_write(this->raw, filename.c_str());
  }

  /**
   * @brief Set a metadata entry.
   * @param key key
   * @param value value
   */
  void set_key_value(const std::string &key, const std::string &value) const {
    pineappl_grid_set_key_value(this->raw, key.c_str(), value.c_str());
  }

  /**
   * @brief Get a metadata entry.
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

  /**
   * @brief Scale grid with a number.
   * This multiplies all subgrids with the given number.
   * @param s factor
   */
  void scale(const double s) const { pineappl_grid_scale(this->raw, s); }

  /**
   * @brief Optimizes the grid representation for space efficiency.
   */
  void optimize() const { pineappl_grid_optimize(this->raw); }
};

/** @brief The central grid object. */
struct Grid : public GridBase {

  /**
   * @brief Constructor.
   * @param lumi luminosity
   * @param orders orders
   * @param bin_limits bin limits
   * @param key_val additional informations
   */
  Grid(const Lumi &lumi, const std::vector<Order> &orders,
       const std::vector<double> &bin_limits, const KeyVal &key_val)
      : GridBase(nullptr) {
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

  /**
   * @brief Number of orders.
   * @return number of orders
   */
  std::size_t order_count() const {
    return pineappl_grid_order_count(this->raw);
  }

  /**
   * @brief Fill grid for the given parameters.
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
   * @brief Perform a convolution of the grid with PDFs.
   * @param pdg_id hadron ID
   * @param pdf PDF
   * @param xi_ren renormalization scale variation
   * @param xi_fac factorization scale variation
   * @param order_mask order mask
   * @param lumi_mask luminosity mask
   * @return prediction for each bin
   */
  std::vector<double> convolve_with_one(
      const std::int32_t pdg_id, LHAPDF::PDF &pdf, const double xi_ren = 1.0,
      const double xi_fac = 1.0, const std::vector<bool> &order_mask = {},
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
    pineappl_grid_convolve_with_one(this->raw, pdg_id, xfx, alphas, &pdf,
                                    raw_order_mask.get(), raw_lumi_mask.get(),
                                    xi_ren, xi_fac, results.data());
    return results;
  }
};

/** @brief The central grid object in the v1 representation. */
struct GridV1 : public GridBase {
  /**
   * @brief Constructor.
   * @param orders orders
   * @param channels channels
   * @param pid_basis pid_basis
   * @param pids pids
   * @param convolutions_types convolution_types
   * @param interp interp
   * @param bin_limits bin_limits
   * @param mu_scales indexed representing the scales
   */
  GridV1(std::vector<OrderV1> &orders, const Channels &channels,
         PidBasis pid_basis, std::vector<int32_t> pids,
         std::vector<ConvType> &convolution_types,
         std::vector<Kinematics> &kinematics, std::vector<InterpTuples> &interp,
         std::vector<double> &bin_limits, std::vector<std::size_t> &mu_scales)
      : GridBase(nullptr) {
    const std::size_t n_orders = orders.size();
    const std::size_t n_bins = bin_limits.size() - 1;
    const std::size_t n_convs = convolution_types.size();

    // Cast the Orders
    std::vector<std::uint8_t> raw_orders;
    for (const OrderV1 &order : orders) {
      raw_orders.push_back(order.alphas);
      raw_orders.push_back(order.alpha);
      raw_orders.push_back(order.logxir);
      raw_orders.push_back(order.logxif);
      raw_orders.push_back(order.logxia);
    }

    this->raw = pineappl_grid_new2(
        pid_basis, channels.raw, n_orders, raw_orders.data(), n_bins,
        bin_limits.data(), n_convs, convolution_types.data(), pids.data(),
        kinematics.data(), interp.data(), mu_scales.data());
    //
  }

  /**
   * @brief Fill grid for the given parameters.
   * @param order perturbative order
   * @param observable value of the observable
   * @param channel index of the channel
   * @param ntuples tuples containing the kinematics {q2, x1, ...}
   * @param weight value of the weighting factor
   */
  void fill(const std::size_t order, const double observable,
            const std::size_t channel, std::vector<double> &ntuples,
            const double weight) const {
    pineappl_grid_fill2(this->raw, order, observable, channel, ntuples.data(),
                        weight);
  }
};

}  // namespace PineAPPL

#endif  // PineAPPL_HPP_
