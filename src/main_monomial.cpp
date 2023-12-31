#include "monomial.hpp"
#include <iostream>

using namespace geometry;

int main()
{
    Monomial2D monomial2D(2, 3, 5.0);
    Monomial2D monomial2D_1(3, 1, 2.0);
    Monomial3D monomial3D(1, 0, 2, 7.5);
    Monomial<4> monomial4D({1, 2, 3, 4}, 5);
    Monomial<2> monomial2D_2({1, 4}, 2.5);

    // Access exponents and coefficient
    std::cout << "Exponent of x of monomial2D is " << monomial2D.getExponents()[0] << std::endl;
    std::cout << "Coefficient of monomial2D is " << monomial2D.getCoefficient() << std::endl;

    std::cout << "Monomial2D: " << monomial2D << std::endl;
    std::cout << "Monomial2D_1: " << monomial2D_1 << std::endl;
    std::cout << "Monomial3D: " << monomial3D << std::endl;
    std::cout << "Monomial4D: " << monomial4D << std::endl;
    std::cout << "Monomial2D_2, base instance: " << monomial2D_2 << std::endl;

    // Product
    auto product = monomial2D * monomial2D_1;
    std::cout << "Product of monomial2D and monomial2D_1: " << product << std::endl;
    auto product3D = monomial3D * monomial3D;
    std::cout << "monomial3D squared: " << product3D << std::endl;
    auto product4D = monomial4D * monomial4D;
    std::cout << "monomial4D squared: " << product4D << std::endl;

    // Derivative
    auto derivative = monomial2D.derivative(1); // Derivative with respect to x2
    std::cout << "Derivative of monomial2D with respect to y: " << derivative << " or called with .dy() " << monomial2D.dy() << std::endl;
    auto derivative3D = monomial3D.derivative(0); // Derivative with respect to x1
    std::cout << "Derivative of monomial3D with respect to x: " << derivative3D << " or called with .dx() " << monomial3D.dx() << std::endl;

    // Compute monomials in 2D up to order k
    monomial2D.computeMonomialsUpToOrder(4);
    std::cout<<"getting monomials in 2D up to order 2"<<std::endl;
    std::vector<Monomial2D> monomialOrdered2D = Monomial2D::getMonomialsOrdered(2);
    for (const auto &m : monomialOrdered2D)
    {
        std::cout << m << std::endl;
    }
    std::cout << monomialOrdered2D.size() << std::endl;

    // Compute monomials in 3D up to order k
    monomial3D.computeMonomialsUpToOrder(4);
    std::cout<<"getting monomials in 3D up to order 7"<<std::endl;
    std::vector<Monomial3D> monomialOrdered3D = Monomial3D::getMonomialsOrdered(7);
    for (const auto &m : monomialOrdered3D)
    {
        std::cout << m << std::endl;
    }
    std::cout << monomialOrdered3D.size() << std::endl;

    // Evaluate monomial
    std::cout << monomial3D.evaluate(Point3D(1, 2, 3)) << std::endl;

    // Get order of a monomial
    std::cout << monomial2D.getOrder() << std::endl;
    std::cout << monomial3D.getOrder() << std::endl;
    std::cout << monomial4D.getOrder() << std::endl;

    // Test polynomial
    std::cout << "polynomial 1:" << std::endl;
    Polynomial<3> polynomial1({Monomial3D(1, 1, 1, 1.0), Monomial3D(2, 2, 2, 1.0), Monomial3D(2, 2, 2, 3.0)});
    for (const auto &monomialPair : polynomial1.getPolynomial())
    {
        std::cout << monomialPair.second << std::endl;
    }
    std::cout << "polynomial 2:" << std::endl;
    Polynomial<3> polynomial2({Monomial3D(1, 2, 3, 1.0)});
    for (const auto &monomialPair : polynomial2.getPolynomial())
    {
        std::cout << monomialPair.second << std::endl;
    }
    std::cout << "product:" << std::endl;
    auto polyprod = polynomial1 * polynomial2;
    for (const auto &monomialPair : polyprod.getPolynomial())
    {
        std::cout << monomialPair.second << std::endl;
    }
    std::cout<<std::endl;

    // Test trinomial power
    std::cout<<"Trinomial power test"<<std::endl;
    LinearTrinomialPower tp(1, 1, 1, 3);
    for (const auto &monomialPair : tp.getPolynomial())
    {
        std::cout << monomialPair.second << std::endl;
    }
    std::cout<<std::endl;

    // Test toPolynomial2D
    std::cout<<"toPolynomial2D test"<<std::endl;
    Monomial3D m3D(0, 0, 2, 1.0);
    auto m3Din2D = toPolynomial2D(m3D, sqrt(3)/2, sqrt(2)/2, Point3D(0.25, 0.25, 0.25), Point3D(0.5, 0.25, 0.25), Point3D(0, 0, 1), Point3D(0, 1, 0), Point3D(-1, 0, 0));
    for (const auto &monomialPair : m3Din2D.getPolynomial())
    {
        std::cout << monomialPair.second << std::endl;
    }
    std::cout<<"order of m3D = "<<m3D.getOrder()<<", order of m3Din2D = "<<m3Din2D.getOrder()<<std::endl;

    /*
    // Test getLaplaciansToMonomialsOrdered in 2D
    auto laplacians2D=Monomial2D::getLaplaciansToMonomialsOrdered(2);
    for (const auto& l:laplacians2D)
    {
        std::cout<<l.first.first<<" "<<monomialOrdered2D[l.first.second]<<std::endl;
        std::cout<<l.second.first<<" "<<monomialOrdered2D[l.second.second]<<std::endl;
    }

    monomialOrdered2D=Monomial2D::getMonomialsOrdered(5);
    laplacians2D=Monomial2D::getLaplaciansToMonomialsOrdered(5);
    for (const auto& l:laplacians2D)
    {
        std::cout<<l.first.first<<" "<<monomialOrdered2D[l.first.second]<<std::endl;
        std::cout<<l.second.first<<" "<<monomialOrdered2D[l.second.second]<<std::endl;
    }

    // Test getLaplaciansToMonomialsOrdered in 3D
    std::cout<<"Test getLaplaciansToMonomialsOrdered in 3D"<<std::endl;
    auto laplacians3D=Monomial3D::getLaplaciansToMonomialsOrdered(3);
    for (const auto& l:laplacians3D)
    {
        std::cout<<l[0].first<<" "<<monomialOrdered3D[l[0].second]<<std::endl;
        std::cout<<l[1].first<<" "<<monomialOrdered3D[l[1].second]<<std::endl;
        std::cout<<l[2].first<<" "<<monomialOrdered3D[l[2].second]<<std::endl;
        std::cout<<std::endl;
    }

    // Test getLaplaciansToMonomialsOrdered in 3D
    std::cout<<"Test getGradientsToMonomialsOrdered in 3D"<<std::endl;
    auto gradients3D=Monomial3D::getGradientsToMonomialsOrdered(3);
    for (const auto& g:gradients3D)
    {
        std::cout<<g[0].first<<" "<<monomialOrdered3D[g[0].second]<<std::endl;
        std::cout<<g[1].first<<" "<<monomialOrdered3D[g[1].second]<<std::endl;
        std::cout<<g[2].first<<" "<<monomialOrdered3D[g[2].second]<<std::endl;
        std::cout<<std::endl;
    }
*/
    return 0;
}