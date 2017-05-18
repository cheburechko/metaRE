#include <algorithm>
#include <stdexcept>
#include <boost/math/distributions/hypergeometric.hpp>

enum class Alternative{LESS, GREATER, TWO_SIDED};

Alternative strToAlternative(std::string val) {
    if (!val.compare("two.sided")) {
        return Alternative::TWO_SIDED;
    } else if (!val.compare("less")) {
        return Alternative::LESS;
    } else if (!val.compare("greater")) {
        return Alternative::GREATER;
    } else {
        throw std::invalid_argument("Illegal alternative");
    }
}

int find_k(int left, int right, int sign, double value, boost::math::hypergeometric& hgd) {
    if (right - left < 2) {
        double leftval = boost::math::pdf(hgd, left);
        double rightval = boost::math::pdf(hgd, right);
        if (leftval <= value && rightval <= value) {
            if (leftval > rightval) {
                return left;
            } else if (leftval < rightval){
                return right;
            } else {
                return sign == 1 ? left : right;
            }
        } else if (leftval <= value) {
            return left;
        } else {
            return right;
        }
    }
    int mid = (left+right)/2;
    double mid_value = boost::math::pdf(hgd, mid);
    if (value == mid_value) {
        return mid;
    }
    if (sign * (value - mid_value) > 0) {
        return find_k(mid, right, sign, value, hgd);
    } else {
        return find_k(left, mid, sign, value, hgd);
    }
}

double fisherTest(double eff1, double n1, double eff2, double n2, Alternative alternative){
    int N = n1 + n2;
    int eff = eff1 + eff2;

    boost::math::hypergeometric hgd(eff, n1, N);

    if (alternative ==  Alternative::TWO_SIDED) {
        int mode = double(eff+1) * (n1+1) / (N+2);

        int min_eff = std::max(0., eff-n2);
        int max_eff = std::min(eff, int(n1));

        int min_k, max_k;
        double cutoff = boost::math::pdf(hgd, eff1);
        if (eff1 < mode) {
            if (boost::math::pdf(hgd, max_eff) > cutoff) {
                return boost::math::cdf(hgd, eff1);
            }
            min_k = eff1;
            max_k = find_k(mode, max_eff, -1, cutoff, hgd);
        } else {
            if (boost::math::pdf(hgd, min_eff) > cutoff) {
                return boost::math::cdf(complement(hgd, eff1-1));
            }
            max_k = eff1;
            min_k = find_k(min_eff, mode, 1, cutoff, hgd);
        }
        if (max_k == min_k) max_k++;

        double p = boost::math::cdf(hgd, min_k) + boost::math::cdf(complement(hgd, max_k-1));

        return p;
    } else if (alternative == Alternative::GREATER) {
        if (int(eff1) == int(eff - n2) || eff1 == 0) {
            return 1.;
        }
        return boost::math::cdf(complement(hgd, eff1-1));
    } else if (alternative == Alternative::LESS) {
        return boost::math::cdf(hgd, eff1);
    } else {
        throw std::invalid_argument("Illegal alternative");
    }
}
