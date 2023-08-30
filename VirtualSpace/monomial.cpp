#include "monomial.hpp"

unsigned int Monomial2D::order = 0;
std::vector<Monomial2D> Monomial2D::monomials_ordered = {Monomial2D(0, 0, 1.0)};
std::vector<std::pair<std::pair<real, std::size_t>,
                      std::pair<real, std::size_t>>>
    Monomial2D::laplacians_to_monomials_ordered = {};

void Monomial2D::computeMonomialsUpToOrder(unsigned int order_)
{
    if (order_ > order)
    {
        std::size_t j;
        for (std::size_t kk = (order + 1); kk <= order_; kk++)
        {
            for (std::size_t i = kk; i > 0; i--)
            {
                j = kk - i;
                monomials_ordered.emplace_back(i, j, 1.0);
            }
            monomials_ordered.emplace_back(0, kk, 1.0);
        }
        order = order_;
    }
};

const std::vector<Monomial2D> Monomial2D::getMonomialsOrdered(unsigned int order_)
{
    if (order_ > order) // Check if the cache for n doesn't exist
    {
        computeMonomialsUpToOrder(order_); // Initialize the cache for n
    }
    std::size_t endIndex = (order_ + 1) * (order_ + 2) / 2;
    return std::vector<Monomial2D>(monomials_ordered.begin(), monomials_ordered.begin() + endIndex);
}

void Monomial2D::computeLaplaciansToMonomialsOrdered(unsigned int order_)
{
    auto monomials = getMonomialsOrdered(order_);
    std::size_t n_kminus2 = 0;
    if (order_ > 1)
        n_kminus2 = getMonomialsOrdered(order_ - 2).size();
    for (std::size_t i = laplacians_to_monomials_ordered.size(); i < monomials.size(); i++)
    {
        auto mdxx = monomials[i].dx().dx();
        auto mdyy = monomials[i].dy().dy();
        real c1(0.0), c2(0.0);
        std::size_t m1(0), m2(0);
        if (mdxx.getCoefficient() != 0.0)
        {
            for (std::size_t j = 0; j < n_kminus2; j++)
            {
                if (monomials[j].getExponents() == mdxx.getExponents())
                {
                    c1 = mdxx.getCoefficient();
                    m1 = j;
                    break;
                }
            }
        }
        if (mdyy.getCoefficient() != 0.0)
        {
            for (std::size_t j = 0; j < n_kminus2; j++)
            {
                if (monomials[j].getExponents() == mdyy.getExponents())
                {
                    c2 = mdyy.getCoefficient();
                    m2 = j;
                    break;
                }
            }
        }
        laplacians_to_monomials_ordered.emplace_back(std::make_pair(std::make_pair(c1, m1),
                                                                    std::make_pair(c2, m2)));
    }
}

const std::vector<std::pair<std::pair<real, std::size_t>,
                            std::pair<real, std::size_t>>>
Monomial2D::getLaplaciansToMonomialsOrdered(unsigned int order_)
{
    if ((laplacians_to_monomials_ordered.size() < ((order_ + 1) * (order_ + 2) / 2))) // Check if the cache for n doesn't exist
    {
        computeLaplaciansToMonomialsOrdered(order_); // Initialize the cache for n
    }
    std::size_t endIndex = (order_ + 1) * (order_ + 2) / 2;
    return std::vector<std::pair<std::pair<real, std::size_t>,
                                 std::pair<real, std::size_t>>>(laplacians_to_monomials_ordered.begin(), laplacians_to_monomials_ordered.begin() + endIndex);
}


unsigned int Monomial3D::order = 0;
std::vector<Monomial3D> Monomial3D::monomials_ordered = {Monomial3D(0, 0, 0, 1.0)};
std::vector<std::array<std::pair<real, std::size_t>, 3>>
    Monomial3D::gradients_to_monomials_ordered = {};
std::vector<std::array<std::pair<real, std::size_t>, 3>>
    Monomial3D::laplacians_to_monomials_ordered = {};

void Monomial3D::computeMonomialsUpToOrder(unsigned int order_)
{
    if (order_ > order)
    {
        std::size_t j;
        for (std::size_t kk = (order + 1); kk <= order_; kk++)
        {
            for (std::size_t k = 0; k < (kk / 3 + 1); k++)
            {
                for (std::size_t i = (kk - 2 * k); i > k; i--)
                {
                    j = kk - i - k;
                    monomials_ordered.emplace_back(i, j, k, 1.0);
                }
                for (std::size_t i = (kk - 2 * k); i > k; i--)
                {
                    j = kk - i - k;
                    monomials_ordered.emplace_back(k, i, j, 1.0);
                }
                for (std::size_t i = (kk - 2 * k); i > k; i--)
                {
                    j = kk - i - k;
                    monomials_ordered.emplace_back(j, k, i, 1.0);
                }
            }
            if (kk % 3 == 0)
            {
                monomials_ordered.emplace_back(kk / 3, kk / 3, kk / 3, 1.0);
            }
        }
    }
    order = order_;
}

const std::vector<Monomial3D>
Monomial3D::getMonomialsOrdered(unsigned int order_)
{
    if (order_ > order) // Check if the cache for n doesn't exist
    {
        computeMonomialsUpToOrder(order_); // Initialize the cache for n
    }
    std::size_t endIndex = (order_ + 1) * (order_ + 2) * (order_ + 3) / 6;
    return std::vector<Monomial3D>(monomials_ordered.begin(), monomials_ordered.begin() + endIndex);
}

void Monomial3D::computeGradientsToMonomialsOrdered(unsigned int order_)
{
    auto monomials = getMonomialsOrdered(order_);
    std::size_t nu_kminus1 = 0;
    if (order_ > 0)
        nu_kminus1 = getMonomialsOrdered(order_ - 1).size();
    for (std::size_t i = gradients_to_monomials_ordered.size(); i < monomials.size(); i++)
    {
        auto mdx = monomials[i].dx();
        auto mdy = monomials[i].dy();
        auto mdz = monomials[i].dz();
        real c1(0.0), c2(0.0), c3(0.0);
        std::size_t m1(0), m2(0), m3(0);
        if (mdx.getCoefficient() != 0.0)
        {
            for (std::size_t j = 0; j < nu_kminus1; j++)
            {
                if (monomials[j].getExponents() == mdx.getExponents())
                {
                    c1 = mdx.getCoefficient();
                    m1 = j;
                    break;
                }
            }
        }
        if (mdy.getCoefficient() != 0.0)
        {
            for (std::size_t j = 0; j < nu_kminus1; j++)
            {
                if (monomials[j].getExponents() == mdy.getExponents())
                {
                    c2 = mdy.getCoefficient();
                    m2 = j;
                    break;
                }
            }
        }
        if (mdz.getCoefficient() != 0.0)
        {
            for (std::size_t j = 0; j < nu_kminus1; j++)
            {
                if (monomials[j].getExponents() == mdz.getExponents())
                {
                    c3 = mdz.getCoefficient();
                    m3 = j;
                    break;
                }
            }
        }
        gradients_to_monomials_ordered.emplace_back(std::array<std::pair<real, std::size_t>, 3>{std::make_pair(c1, m1),
                                                                                                std::make_pair(c2, m2),
                                                                                                std::make_pair(c3, m3)});
    }
}

const std::vector<std::array<std::pair<real, std::size_t>, 3>>
Monomial3D::getGradientsToMonomialsOrdered(unsigned int order_)
{
    if ((gradients_to_monomials_ordered.size() < ((order_ + 1) * (order_ + 2) + (order_ + 3) / 6))) // Check if the cache for n doesn't exist
    {
        computeGradientsToMonomialsOrdered(order_); // Initialize the cache for n
    }
    std::size_t endIndex = (order_ + 1) * (order_ + 2) * (order_ + 3) / 6;
    return std::vector<std::array<std::pair<real, std::size_t>, 3>>(gradients_to_monomials_ordered.begin(), gradients_to_monomials_ordered.begin() + endIndex);
}

void Monomial3D::computeLaplaciansToMonomialsOrdered(unsigned int order_)
{
    auto monomials = getMonomialsOrdered(order_);
    std::size_t nu_kminus2 = 0;
    if (order_ > 1)
        nu_kminus2 = getMonomialsOrdered(order_ - 2).size();
    for (std::size_t i = laplacians_to_monomials_ordered.size(); i < monomials.size(); i++)
    {
        auto mdxx = monomials[i].dx().dx();
        auto mdyy = monomials[i].dy().dy();
        auto mdzz = monomials[i].dz().dz();
        real c1(0.0), c2(0.0), c3(0.0);
        std::size_t m1(0), m2(0), m3(0);
        if (mdxx.getCoefficient() != 0.0)
        {
            for (std::size_t j = 0; j < nu_kminus2; j++)
            {
                if (monomials[j].getExponents() == mdxx.getExponents())
                {
                    c1 = mdxx.getCoefficient();
                    m1 = j;
                    break;
                }
            }
        }
        if (mdyy.getCoefficient() != 0.0)
        {
            for (std::size_t j = 0; j < nu_kminus2; j++)
            {
                if (monomials[j].getExponents() == mdyy.getExponents())
                {
                    c2 = mdyy.getCoefficient();
                    m2 = j;
                    break;
                }
            }
        }
        if (mdzz.getCoefficient() != 0.0)
        {
            for (std::size_t j = 0; j < nu_kminus2; j++)
            {
                if (monomials[j].getExponents() == mdzz.getExponents())
                {
                    c3 = mdzz.getCoefficient();
                    m3 = j;
                    break;
                }
            }
        }
        laplacians_to_monomials_ordered.emplace_back(std::array<std::pair<real, std::size_t>, 3>{std::make_pair(c1, m1),
                                                                                                 std::make_pair(c2, m2),
                                                                                                 std::make_pair(c3, m3)});
    }
}

const std::vector<std::array<std::pair<real, std::size_t>, 3>>
Monomial3D::getLaplaciansToMonomialsOrdered(unsigned int order_)
{
    if ((laplacians_to_monomials_ordered.size() < ((order_ + 1) * (order_ + 2) + (order_ + 3) / 6))) // Check if the cache for n doesn't exist
    {
        computeLaplaciansToMonomialsOrdered(order_); // Initialize the cache for n
    }
    std::size_t endIndex = (order_ + 1) * (order_ + 2) * (order_ + 3) / 6;
    return std::vector<std::array<std::pair<real, std::size_t>, 3>>(laplacians_to_monomials_ordered.begin(), laplacians_to_monomials_ordered.begin() + endIndex);
}
