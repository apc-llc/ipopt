// Copyright (C) 2004, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: IpStdCInterface.cpp 2400 2013-10-19 18:38:36Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpStdCInterface.h"
#include "IpStdInterfaceTNLP.hpp"
#include "IpOptionsList.hpp"
#include "IpIpoptApplication.hpp"

#include <vector>

using namespace Ipopt;
using namespace std;

struct IpoptProblemInfo
{
  Index n;
  vector<Number> x_L;
  vector<Number> x_U;
  Index m;
  vector<Number> g_L;
  vector<Number> g_U;
  Index nele_jac;
  Index nele_hess;
  Index index_style;
  Eval_F_CB eval_f;
  Eval_G_CB eval_g;
  Eval_Grad_F_CB eval_grad_f;
  Eval_Jac_G_CB eval_jac_g;
  Eval_H_CB eval_h;
  Intermediate_CB intermediate_cb;
  SmartPtr<IpoptApplication> app;
  Number obj_scaling;
  vector<Number> x_scaling;
  vector<Number> g_scaling;

  SmartPtr<TNLP> tnlp;
  vector<Number> start_x;
  vector<Number> start_lam;
  vector<Number> start_z_L;
  vector<Number> start_z_U;

  IpoptProblemInfo(
    Index n_,
    Number* x_L_,
    Number* x_U_,
    Index m_,
    Number* g_L_,
    Number* g_U_,
    Index nele_jac_,
    Index nele_hess_,
    Index index_style_,
    Eval_F_CB eval_f_,
    Eval_G_CB eval_g_,
    Eval_Grad_F_CB eval_grad_f_,
    Eval_Jac_G_CB eval_jac_g_,
    Eval_H_CB eval_h_) :

	n(n_),
	x_L(n),
	x_U(n),
	m(m_),
	g_L(m),
	g_U(m),
	nele_jac(nele_jac_),
	nele_hess(nele_hess_),
	index_style(index_style_),
	eval_f(eval_f_),
	eval_g(eval_g_),
	eval_grad_f(eval_grad_f_),
	eval_jac_g(eval_jac_g_),
	eval_h(eval_h_),
	intermediate_cb(NULL),
	app(new IpoptApplication()),
	obj_scaling(1),
	x_scaling(n),
	g_scaling(m)

  {
    for (Index i=0; i<n; i++) {
      x_L[i] = x_L_[i];
    }
    for (Index i=0; i<n; i++) {
      x_U[i] = x_U_[i];
    }

    if (m>0) {
      for (Index i=0; i<m; i++) {
        g_L[i] = g_L_[i];
      }
      for (Index i=0; i<m; i++) {
        g_U[i] = g_U_[i];
      }
    }
  }
};

IpoptProblem CreateIpoptProblem(
  Index n,
  Number* x_L,
  Number* x_U,
  Index m,
  Number* g_L,
  Number* g_U,
  Index nele_jac,
  Index nele_hess,
  Index index_style,
  Eval_F_CB eval_f,
  Eval_G_CB eval_g,
  Eval_Grad_F_CB eval_grad_f,
  Eval_Jac_G_CB eval_jac_g,
  Eval_H_CB eval_h)
{
  // make sure input is Ok
  if (n<1 || m<0 || !x_L || !x_U || (m>0 && (!g_L || !g_U)) ||
      (m==0 && nele_jac != 0) || (m>0 && nele_jac < 1) || nele_hess < 0 ||
      !eval_f || !eval_grad_f || (m>0 && (!eval_g || !eval_jac_g))) {
    return NULL;
  }

  IpoptProblem retval = new IpoptProblemInfo(
    n,
    x_L,
    x_U,
    m,
    g_L,
    g_U,
    nele_jac,
    nele_hess,
    index_style,
    eval_f,
    eval_g,
    eval_grad_f,
    eval_jac_g,
    eval_h);

  retval->app->RethrowNonIpoptException(false);

  return retval;
}

void FreeIpoptProblem(IpoptProblem ipopt_problem)
{
  ipopt_problem->app = NULL;

  delete ipopt_problem;
}

Bool AddIpoptStrOption(IpoptProblem ipopt_problem, char* keyword, char* val)
{
  string tag(keyword);
  string value(val);
  return (Bool) ipopt_problem->app->Options()->SetStringValue(tag, value);
}

Bool AddIpoptNumOption(IpoptProblem ipopt_problem, char* keyword, Number val)
{
  string tag(keyword);
  Number value=val;
  return (Bool) ipopt_problem->app->Options()->SetNumericValue(tag, value);
}

Bool AddIpoptIntOption(IpoptProblem ipopt_problem, char* keyword, Int val)
{
  string tag(keyword);
  Index value=val;
  return (Bool) ipopt_problem->app->Options()->SetIntegerValue(tag, value);
}

Bool OpenIpoptOutputFile(IpoptProblem ipopt_problem, char* file_name,
                         Int print_level)
{
  string name(file_name);
  EJournalLevel level = EJournalLevel(print_level);
  return (Bool) ipopt_problem->app->OpenOutputFile(name, level);
}

Bool SetIpoptProblemScaling(IpoptProblem ipopt_problem,
                            Number obj_scaling,
                            Number* x_scaling,
                            Number* g_scaling)
{
  ipopt_problem->obj_scaling = obj_scaling;
  if (x_scaling) {
    for (Index i=0; i<ipopt_problem->n; i++) {
      ipopt_problem->x_scaling[i] = x_scaling[i];
    }
  }
  if (g_scaling) {
    for (Index i=0; i<ipopt_problem->m; i++) {
      ipopt_problem->g_scaling[i] = g_scaling[i];
    }
  }

  return (Bool)true;
}

Bool SetIntermediateCallback(IpoptProblem ipopt_problem,
                             Intermediate_CB intermediate_cb)
{
  ipopt_problem->intermediate_cb = intermediate_cb;
  return (Bool)true;
}

enum ::ApplicationReturnStatus IpoptSolve(
  IpoptProblem ipopt_problem,
  Number* x,
  Number* g,
  Number* obj_val,
  Number* mult_g,
  Number* mult_x_L,
  Number* mult_x_U,
  UserDataPtr user_data)
{
  Ipopt::ApplicationReturnStatus status;

  // Create the original nlp, if not already done.
  if (GetRawPtr(ipopt_problem->tnlp) == NULL) {
    // Initialize and process options
    Ipopt::ApplicationReturnStatus retval = ipopt_problem->app->Initialize();
    if (retval!=Ipopt::Solve_Succeeded) {
      return (::ApplicationReturnStatus) retval;
    }

    if (!x) {
      ipopt_problem->app->Jnlst()->Printf(J_ERROR, J_MAIN,
                                          "Error: Array x with starting point information is NULL.");
      return (::ApplicationReturnStatus) ::Invalid_Problem_Definition;
    }

    // Copy the starting point information
    ipopt_problem->start_x.resize(ipopt_problem->n);
    for (Index i=0; i<ipopt_problem->n; i++) {
      ipopt_problem->start_x[i] = x[i];
    }
    if (mult_g) {
      ipopt_problem->start_lam.resize(ipopt_problem->m);
      for (Index i=0; i<ipopt_problem->m; i++) {
        ipopt_problem->start_lam[i] = mult_g[i];
      }
    }
    if (mult_x_L) {
      ipopt_problem->start_z_L.resize(ipopt_problem->n);
      for (Index i=0; i<ipopt_problem->n; i++) {
        ipopt_problem->start_z_L[i] = mult_x_L[i];
      }
    }
    if (mult_x_U) {
      ipopt_problem->start_z_U.resize(ipopt_problem->n);
      for (Index i=0; i<ipopt_problem->n; i++) {
        ipopt_problem->start_z_U[i] = mult_x_U[i];
      }
    }

    try {
      ipopt_problem->tnlp = new StdInterfaceTNLP(ipopt_problem->n,
                                                 &ipopt_problem->x_L[0],
                                                 &ipopt_problem->x_U[0],
                                                 ipopt_problem->m,
                                                 (ipopt_problem->m > 0) ? &ipopt_problem->g_L[0] : NULL,
                                                 (ipopt_problem->m > 0) ? &ipopt_problem->g_U[0] : NULL,
                                                 ipopt_problem->nele_jac,
                                                 ipopt_problem->nele_hess,
                                                 ipopt_problem->index_style,
                                                 &ipopt_problem->start_x[0],
                                                 &ipopt_problem->start_lam[0],
                                                 &ipopt_problem->start_z_L[0],
                                                 &ipopt_problem->start_z_U[0],
                                                 ipopt_problem->eval_f,
                                                 ipopt_problem->eval_g,
                                                 ipopt_problem->eval_grad_f,
                                                 ipopt_problem->eval_jac_g,
                                                 ipopt_problem->eval_h,
                                                 ipopt_problem->intermediate_cb,
                                                 x, mult_x_L, mult_x_U, g, mult_g,
                                                 obj_val, user_data,
                                                 ipopt_problem->obj_scaling,
                                                 &ipopt_problem->x_scaling[0],
                                                 (ipopt_problem->m > 0) ? &ipopt_problem->g_scaling[0] : NULL);
    }
    catch (INVALID_STDINTERFACE_NLP& exc) {
      exc.ReportException(*ipopt_problem->app->Jnlst(), J_ERROR);
      status = Ipopt::Invalid_Problem_Definition;
    }
  }

  try {
    status = ipopt_problem->app->OptimizeTNLP(ipopt_problem->tnlp);
  }
  catch( IpoptException& exc ) {
    exc.ReportException(*ipopt_problem->app->Jnlst(), J_ERROR);
    status = Ipopt::Unrecoverable_Exception;
  }

  return (::ApplicationReturnStatus) status;
}

