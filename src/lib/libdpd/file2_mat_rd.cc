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

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include <libpsio/psio.h>
#include "dpd.h"

namespace psi {

int DPD::file2_mat_rd(dpdfile2 *File)
{
    int h, my_irrep, rowtot, coltot;
    psio_address irrep_ptr, next_address;

    my_irrep = File->my_irrep;

    if(File->incore) return 0; /* We already have this data in core */

    /* If data doesn't actually exist on disk, we just leave */
    if(psio_->tocscan(File->filenum, File->label) == NULL) return 1;

    for(h=0; h < File->params->nirreps; h++) {
        irrep_ptr = File->lfiles[h];
        rowtot = File->params->rowtot[h];
        coltot = File->params->coltot[h^my_irrep];

        if(rowtot && coltot)
            psio_->read(File->filenum, File->label, (char *) File->matrix[h][0],
                    rowtot*coltot*sizeof(double), irrep_ptr, &next_address);
    }

    return 0;
}

}
