/**
 * \file       common.hpp
 * \author     Bryan Flynt
 * \date       Jun 09, 2022
 * \copyright  Copyright (C) 2021 Bryan Flynt - All Rights Reserved
 */
#pragma once

#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

namespace common {

/// Print Coordinates to Console
/**
 */
template <typename T>
void dump_line(const std::vector<T> xyz) {
    for (std::size_t n = 0; n < xyz.size(); ++n) {
        std::cout << std::setw(15) << std::setprecision(8) << std::scientific << xyz[n] << " ";
    }
    std::cout << std::endl;
}

/// Print Coordinates + Variables to Console
/**
 */
template <typename T>
void dump_line(const std::vector<T> xyz, const std::vector<T> val) {
    for (std::size_t n = 0; n < xyz.size(); ++n) {
        std::cout << std::setw(15) << std::setprecision(8) << std::scientific << xyz[n] << " ";
    }
    for (std::size_t n = 0; n < val.size(); ++n) {
        std::cout << std::setw(15) << std::setprecision(8) << std::scientific << val[n] << " ";
    }
    std::cout << std::endl;
}

/// Print Coordinates + Variables to Console
/**
 */
template <typename T, std::size_t N>
void dump_line(const std::vector<T> xyz, const std::size_t nval, const std::array<T,N> val) {
    for (std::size_t n = 0; n < xyz.size(); ++n) {
        std::cout << std::setw(15) << std::setprecision(8) << std::scientific << xyz[n] << " ";
    }
    for (std::size_t n = 0; n < nval; ++n) {
        std::cout << std::setw(15) << std::setprecision(8) << std::scientific << val[n] << " ";
    }
    std::cout << std::endl;
}

/// Fill Vector with Random Numbers
/**
 */
template <typename T, typename A>
void fill_random(const std::pair<T, T> xyz_range, std::vector<T, A>& xyz) {
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<T> unif(std::get<0>(xyz_range), std::get<1>(xyz_range));
    for (auto& ref : xyz) {
        ref = unif(re);
    }
}

///
/**
 * Fill Vector with Linear Function of First Argument
 */
template <typename T, typename A>
void fill_variables(const std::vector<T, A>& xyz, std::vector<T, A>& var) {
    for (std::size_t i = 0; i < var.size(); ++i) {
        var[i] = i + 3.1415;
        for (std::size_t j = 0; j < xyz.size(); j = j + 2) {
            var[i] += 0.13 * xyz[j];
        }
        for (std::size_t j = 1; j < xyz.size(); j = j + 2) {
            var[i] -= 0.10 * xyz[j];
        }
    }
}

}  // namespace common

/// Class to act as std::arrayu with simple math
/**
 * This class is a simple wrapper around std::array
 * with inefficient math functions added.  It would have
 * been easier to add free math functions to std::array
 * but the GCAL code refuses to use them and complains.
 */
template <typename T, std::size_t N>
struct Array {
    using self_type       = Array<T, N>;
    using container_type  = std::array<T, N>;
    using size_type       = typename container_type::size_type;
    using value_type      = typename container_type::value_type;
    using reference       = typename container_type::reference;
    using const_reference = typename container_type::const_reference;
    using iterator        = typename container_type::iterator;
    using const_iterator  = typename container_type::const_iterator;

    constexpr reference operator[](const size_type index) noexcept { return array_[index]; }

    constexpr const_reference operator[](const size_type index) const noexcept { return array_[index]; }

    constexpr self_type& operator+=(const_reference value) noexcept {
        for (size_type i = 0; i < N; ++i) {
            array_[i] += value;
        }
        return *this;
    }

    constexpr self_type& operator-=(const_reference value) noexcept {
        for (size_type i = 0; i < N; ++i) {
            array_[i] -= value;
        }
        return *this;
    }

    constexpr self_type& operator*=(const_reference value) noexcept {
        for (size_type i = 0; i < N; ++i) {
            array_[i] *= value;
        }
        return *this;
    }

    constexpr self_type& operator/=(const_reference value) noexcept {
        for (size_type i = 0; i < N; ++i) {
            array_[i] /= value;
        }
        return *this;
    }

    constexpr self_type& operator+=(self_type const& rhs) noexcept {
        for (size_type i = 0; i < N; ++i) {
            array_[i] += rhs[i];
        }
        return *this;
    }

    constexpr self_type& operator-=(self_type const& rhs) noexcept {
        for (size_type i = 0; i < N; ++i) {
            array_[i] -= rhs[i];
        }
        return *this;
    }

    constexpr self_type& operator*=(self_type const& rhs) noexcept {
        for (size_type i = 0; i < N; ++i) {
            array_[i] *= rhs[i];
        }
        return *this;
    }

    constexpr self_type& operator/=(self_type const& rhs) noexcept {
        for (size_type i = 0; i < N; ++i) {
            array_[i] /= rhs[i];
        }
        return *this;
    }

    constexpr const_iterator cbegin() const noexcept {
        return array_.cbegin();
    }

    constexpr const_iterator cend() const noexcept {
        return array_.cend();
    }

    constexpr size_type size() const noexcept { return array_.size(); }

    constexpr void fill(value_type const& v) noexcept { array_.fill(v); }

    container_type array_;
};

template <typename T, std::size_t N>
Array<T, N> operator+(const Array<T, N>& lhs, const Array<T, N>& rhs) {
    Array<T, N> ans(lhs);
    ans += rhs;
    return ans;
}

template <typename T, std::size_t N>
Array<T, N> operator*(typename Array<T, N>::value_type const& lhs, const Array<T, N>& rhs) {
    Array<T, N> ans(rhs);
    ans *= lhs;
    return ans;
}

template <typename T, std::size_t N>
Array<T, N> operator*(const Array<T, N>& lhs, typename Array<T, N>::value_type const& rhs) {
    Array<T, N> ans(lhs);
    ans *= rhs;
    return ans;
}

template <typename T, std::size_t N>
Array<T, N> operator/(typename Array<T, N>::value_type const& lhs, const Array<T, N>& rhs) {
    Array<T, N> ans;
    ans.fill(lhs);
    ans /= rhs;
    return ans;
}

template <typename T, std::size_t N>
Array<T, N> operator/(const Array<T, N>& lhs, typename Array<T, N>::value_type const& rhs) {
    Array<T, N> ans(lhs);
    ans /= rhs;
    return ans;
}
