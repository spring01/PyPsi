/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*!
  \file
  \ingroup CHKPT
*/

#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int Chkpt::rd_iopen(void)
{
        int iopen;
        char *keyword;
        keyword = build_keyword("Iopen");

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &iopen, sizeof(int));

        free(keyword);
        return iopen;
}

void Chkpt::wt_iopen(int iopen)
{
        char *keyword;
        keyword = build_keyword("Iopen");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) &iopen, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_iopen()
** Reads in dimensionality of ALPHA and BETA vectors of two-electron
** coupling coefficients for open shells.
**
** Note : IOPEN = MM * (MM + 1), where MM is the total number
** of irreps containing singly occupied orbitals.
**
** returns:
**   iopen = dimensionality of ALPHA and BETA vectors of coupling
**           coefficients for open shells.
** \ingroup CHKPT
*/
        int chkpt_rd_iopen(void)
        {
                return _default_chkpt_lib_->rd_iopen();
        }


/*!
** void chkpt_wt_iopen(int)
** Writes out the dimensionality of ALPHA and BETA vectors of two-electron
** coupling coefficients for open shells.
**
** Note : IOPEN = MM * (MM + 1), where MM is the total number
** of irreps containing singly occupied orbitals.
**
**  arguments:
**   \param iopen = dimensionality of ALPHA and BETA vectors of coupling
**                  coefficients for open shells.
** \ingroup CHKPT
*/
        void chkpt_wt_iopen(int iopen)
        {
                _default_chkpt_lib_->wt_iopen(iopen);
        }
}

