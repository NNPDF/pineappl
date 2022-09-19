/**
 * @file PineAPPL.hpp
 * @brief object oriented interface to PineAPPL
 */
#ifndef PineAPPL_HPP_
#define PineAPPL_HPP_

#include <string>
#include <list>
#include <LHAPDF/LHAPDF.h>

#include <pineappl_capi.h>

/** @brief object oriented interface to PineAPPL  */
namespace PineAPPL {

/** @brief Key-value storage for passing optional information during grid creation */
struct KeyVal {
    /** @brief underlying raw object */ 
    pineappl_keyval *raw;

    /** @brief constructor */
    KeyVal() { this->raw = pineappl_keyval_new(); }

    /** @brief destructor */
    ~KeyVal() { pineappl_keyval_delete(this->raw); }

    
    /** @name setter */
    ///@{
    void set_double(const std::string& key, const double value) const { pineappl_keyval_set_double(this->raw, key.c_str(), value); }
    void set_bool(const std::string& key, const bool value) const { pineappl_keyval_set_bool(this->raw, key.c_str(), value); }
    void set_int(const std::string& key, const int value) const { pineappl_keyval_set_int(this->raw, key.c_str(), value); }
    void set_string(const std::string& key, const std::string& value) const { pineappl_keyval_set_string(this->raw, key.c_str(), value.c_str()); }
    ///@}

    /** @name getter */
    ///@{
    double get_double(const std::string& key) const { return pineappl_keyval_double(this->raw, key.c_str()); }
    bool get_bool(const std::string& key) const { return pineappl_keyval_bool(this->raw, key.c_str()); }
    int get_int(const std::string& key) const { return pineappl_keyval_int(this->raw, key.c_str()); }
    std::string get_string(const std::string& key) const { return pineappl_keyval_string(this->raw, key.c_str()); }
    ///@}

};

/** @brief Entry in luminosity function  */
struct LumiEntry {
    /** @brief first parton id */
    int32_t pid1;

    /** @brief second parton id */
    int32_t pid2;

    /** @brief relative weight */
    double weight;
};

/** @brief Luminosity function */
struct Lumi {
    /** @brief underlying raw object */
    pineappl_lumi *raw;

    /** @brief constructor */
    Lumi(){ this->raw = pineappl_lumi_new(); }

    /** @brief destructor */
    ~Lumi(){ pineappl_lumi_delete(this->raw); }

    /** @brief number of elements */
    size_t count() const { return pineappl_lumi_count(this->raw); }

    /**
     * @brief add a luminosity function
     * @param c luminosity function
     */
    void add(const std::vector<LumiEntry>& c) const {
        size_t n = 0;
        int32_t *pids = new int32_t;
        double *weights = new double;
        for (auto it = c.cbegin(); it != c.cend(); ++it, ++n) {
            pids[2*n] = it->pid1;
            pids[2*n+1] = it->pid2;
            weights[n] = it->weight;
        }
        pineappl_lumi_add(this->raw, n, pids, weights);
        delete pids;
        delete weights;
    }

    /**
     * @brief Returns the number of combinations of the luminosity function `lumi` for the specified entry.
     * @param entry position in lumi
     */
    size_t combinations(size_t entry) const { return pineappl_lumi_combinations(this->raw, entry); }
};

/** @brief Coupling powers for each grid. */
struct Order {
    /** @brief Exponent of the strong coupling. */
    uint32_t alphas;

    /** @brief Exponent of the electromagnetic coupling. */
    uint32_t alpha;

    /** @brief Exponent of the logarithm of the scale factor of the renomalization scale. */
    uint32_t logxir;

    /** @brief Exponent of the logarithm of the scale factor of the factorization scale. */
    uint32_t logxif;
};

class Grid {
    /** @brief number of orders */
    size_t n_orders = 0;

    /** @brief number of bins */
    size_t n_bins = 0;

    /** @brief number of lumis */
    size_t n_lumis = 0;

public:

    /** @brief underlying raw object */
    pineappl_grid *raw;

    /**
     * @brief constructor
     * @param lumi luminosity
     * @param orders orders
     * @param bins bins
     * @param key_val additional informations
     */
    Grid(const Lumi& lumi, const std::vector<Order>& orders, const std::vector<double>& bins, const KeyVal& key_val) {
        this->n_lumis = lumi.count();
        // cast orders
        uint32_t *raw_orders = new uint32_t[4 * orders.size()];
        for (auto it = orders.cbegin(); it != orders.cend(); ++it, ++this->n_orders){
            raw_orders[4*this->n_orders + 0] = it->alphas;
            raw_orders[4*this->n_orders + 1] = it->alpha;
            raw_orders[4*this->n_orders + 2] = it->logxir;
            raw_orders[4*this->n_orders + 3] = it->logxif;
        }
        //this->n_bins = bins.size();
        // cast bins
        double *raw_bins = new double[bins.size()];
        for (auto it = bins.cbegin(); it != bins.cend(); ++it, ++this->n_bins)
            raw_bins[this->n_bins] = *it;
        // now, we got all we need
        this->raw = pineappl_grid_new(lumi.raw, this->n_orders, raw_orders, this->n_bins - 1, raw_bins, key_val.raw);
        // clean up
        delete[] raw_orders;
        delete[] raw_bins;
    }

    /** @brief destructor */
    ~Grid(){ pineappl_grid_delete(this->raw); }

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
    void fill(double x1, double x2, double q2, size_t order, double observable, size_t lumi, double weight) const {
        pineappl_grid_fill(this->raw, x1, x2, q2, order, observable, lumi, weight);
    }

    /**
     * @brief perform a convolution of the grid with PDFs
     * @param pdg_id particle ID
     * @param pdf PDF
     * @param xi_ren renormalization scale variation
     * @param xi_fac factorization scale variation
     * @return prediction for each bin
     * @return prediction for each bin
     */
    std::vector<double> convolute_with_one(int32_t pdg_id, LHAPDF::PDF *pdf, double xi_ren = 1.0, double xi_fac = 1.0 ) const {
        std::vector<bool> order_mask(this->n_orders, true);
        std::vector<bool> lumi_mask(this->n_lumis, true);
        return this->convolute_with_one(pdg_id, pdf, order_mask, lumi_mask, xi_ren, xi_fac);
    }

    /**
     * @brief perform a convolution of the grid with PDFs
     * @param pdg_id particle ID
     * @param pdf PDF
     * @param order_mask order mask
     * @param lumi_mask luminosity mask
     * @param xi_ren renormalization scale variation
     * @param xi_fac factorization scale variation
     * @return prediction for each bin
     */
    std::vector<double> convolute_with_one(int32_t pdg_id, LHAPDF::PDF *pdf, const std::vector<bool>& order_mask, const std::vector<bool>& lumi_mask, double xi_ren = 1.0, double xi_fac = 1.0) const {
        // prepare LHAPDF stuff
        auto xfx = [](int32_t id, double x, double q2, void* pdf) {
            return static_cast <LHAPDF::PDF*> (pdf)->xfxQ2(id, x, q2);
        };
        auto alphas = [](double q2, void* pdf) {
            return static_cast <LHAPDF::PDF*> (pdf)->alphasQ2(q2);
        };
        // cast order_mask
        bool *raw_order_mask = new bool[this->n_orders];
        {size_t j = 0;
        for (auto it = order_mask.cbegin(); it != order_mask.cend(); ++it, ++j)
            raw_order_mask[j] = *it;}
        bool *raw_lumi_mask = new bool[this->n_lumis];
        {size_t j = 0;
        for (auto it = lumi_mask.cbegin(); it != lumi_mask.cend(); ++it, ++j)
            raw_lumi_mask[j] = *it;}
        std::vector<double> results(this->n_bins - 1);
        pineappl_grid_convolute_with_one(this->raw, pdg_id, xfx, alphas, pdf, raw_order_mask, raw_lumi_mask, xi_ren, xi_fac, results.data());
        delete[] raw_order_mask;
        delete[] raw_lumi_mask;
        return results;
    }

    /**
     * @brief Write grid to file
     * @param filename file name
     */
    void write(std::string filename) const {
        pineappl_grid_write(this->raw, filename.c_str());
    }

    /**
     * @brief Set a metadata entry
     * @param key key
     * @param value value
     */
    void set_key_value(std::string key, std::string value) const {
        pineappl_grid_set_key_value(this->raw, key.c_str(), value.c_str());
    }

    /**
     * @brief Get a metadata entry
     * @param key key
     * @return value
     */
    std::string get_key_value(std::string key) const {
        auto* value = pineappl_grid_key_value(this->raw, key.c_str());
        std::string res(value);
        // delete the allocated object
        pineappl_string_delete(value);
        return res;
    }
};

}


#endif // PineAPPL_HPP_
