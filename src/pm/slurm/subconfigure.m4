[#] start of __file__

AC_DEFUN([PAC_SUBCFG_PREREQ_]PAC_SUBCFG_AUTO_SUFFIX,[
])

AC_DEFUN([PAC_SUBCFG_BODY_]PAC_SUBCFG_AUTO_SUFFIX,[


# the pm_names variable is set by the top level configure
build_slurm=no
build_pm_other=no
for pm_name in $pm_names ; do
    if test "X$pm_name" = "Xslurm" ; then
        build_slurm=yes
    else
        build_pm_other=yes
    fi
done

AS_IF([test "X$build_slurm" = "Xyes"],
    [AC_DEFINE([SLURM_PMI_CLIENT], [1], [Define if using slurm pmi client])])
AM_CONDITIONAL([BUILD_PM_SLURM],[test "X$build_slurm" = "Xyes"])
AM_CONDITIONAL([BUILD_PM_OTHER],[test "X$build_pm_other" = "Xyes"])

AM_COND_IF([BUILD_PM_SLURM],[
AC_MSG_NOTICE([RUNNING CONFIGURE FOR src/pm/slurm])

# check conflict with other pms
AM_COND_IF([BUILD_PM_OTHER], [
    AC_MSG_ERROR([Can not select slurm with other process managers])
])

# sets CPPFLAGS and LDFLAGS
PAC_SET_HEADER_LIB_PATH([slurm])

# with_pmi is set by the top level configure

if test "x$with_pmi" = "xsimple" -o "x$with_pmi" = "xpmi1"; then
    AC_CHECK_HEADER([slurm/pmi.h], [], [AC_MSG_ERROR([could not find slurm/pmi.h.  Configure aborted])])
    AC_CHECK_LIB([pmi], [PMI_Init],
             [PAC_PREPEND_FLAG([-lpmi],[LIBS])
              PAC_PREPEND_FLAG([-lpmi], [WRAPPER_LIBS])],
             [AC_MSG_ERROR([could not find the slurm libpmi.  Configure aborted])])
    AC_CHECK_FUNCS([PMI_Ibarrier PMI_Wait])
elif test "x$with_pmi" = "xpmi2/simple" -o "x$with_pmi" = "xpmi2"; then
    USE_PMI2_API=yes
    AC_CHECK_HEADER([slurm/pmi2.h], [], [AC_MSG_ERROR([could not find slurm/pmi2.h.  Configure aborted])])
    AC_CHECK_LIB([pmi2], [PMI2_Init],
             [PAC_PREPEND_FLAG([-lpmi2],[LIBS])
              PAC_PREPEND_FLAG([-lpmi2], [WRAPPER_LIBS])],
             [AC_MSG_ERROR([could not find the slurm libpmi2.  Configure aborted])])
    AC_CHECK_FUNCS([PMI2_KVS_Ifence PMI2_KVS_Wait])
    AC_CHECK_FUNCS([PMI2_Iallgather PMI2_Iallgather_wait])
    AC_CHECK_FUNCS([PMI2_SHMEM_Iallgather PMI2_SHMEM_Iallgather_wait])
else
    AC_MSG_ERROR([Selected PMI ($with_pmi) is not compatible with slurm])
fi

])dnl end COND_IF

])dnl end BODY macro

[#] end of __file__
