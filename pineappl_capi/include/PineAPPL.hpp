/**
 * @brief Object oriented interface to PineAPPL using only v1 features.
 */
#ifndef PineAPPL_HPP_
#define PineAPPL_HPP_

#include <LHAPDF/LHAPDF.h>
#include <pineappl_capi.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

/** @brief Object oriented interface to PineAPPL.*/
namespace PineAPPL {

  // TODO: Add checks within various functions/calls
  // TODO: Implement `grid_convolve`

/** @brief Entry in sub-channel function. A sub-channel consists of a vector of
 * tuple. Each tuple contains two elements. The first elements store the list of
 * PIDs as a vector while the second represent the weighting factor. Note that
 * the number of PIDs must correspond to the number of convolutions involved.*/
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

/** @brief Extension of `Order` to accommodate for v1 representation. */
struct Order {
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
struct Grid {
  /** @brief Underlying raw object. */
  pineappl_grid *raw;

  /** @brief Constructor (protected to avoid direct instantiation). */
 protected:
  Grid(pineappl_grid *grid) : raw(grid) {}

  /** @brief Deleted copy/move semantics. */
  Grid() = delete;
  Grid(const Grid &) = delete;
  Grid(Grid &&) = delete;

  Grid &operator=(const Grid &) = delete;
  Grid &operator=(Grid &&) = delete;

 public:
  /** @brief Destructor. */
  virtual ~Grid() { pineappl_grid_delete(this->raw); }

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
  Grid(std::vector<Order> &orders, const Channels &channels,
       pineappl_pid_basis pid_basis, std::vector<int32_t> pids,
       std::vector<pineappl_conv_type> &convolution_types,
       std::vector<pineappl_kinematics> &kinematics,
       std::vector<pineappl_interp_tuples> &interp,
       std::vector<double> &bin_limits, std::vector<std::size_t> &mu_scales)
      : Grid(nullptr) {
    const std::size_t n_orders = orders.size();
    const std::size_t n_bins = bin_limits.size() - 1;
    const std::size_t n_convs = convolution_types.size();

    // Cast the Orders
    std::vector<std::uint8_t> raw_orders;
    for (const Order &order : orders) {
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

}  // namespace PineAPPL

#endif  // PineAPPL_HPP_
