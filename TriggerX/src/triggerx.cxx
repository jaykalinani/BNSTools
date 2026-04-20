#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Constants.h"
#include "cctk_Groups.h"
#include "cctk_Parameters.h"
#include "cctk_WarnLevel.h"

#include <AMReX.H>
#include <tuple.hxx>
#include <loop.hxx>

namespace CarpetX {
using Loop::dim;
}

#include <reduction.hxx>

typedef struct {
  int number;
  int debug;
  int *active;
  int *checked_variable;
  int *steered_scalar;
  int *last_checked;
  int *trigger_once;
  int *trigger_count;
  const char **relation;
  const char **reduction;
  const char **checked_parameter_thorn;
  const char **checked_parameter_name;
  const char **reaction;
  CCTK_REAL *checked_value;
} TriggerXGH;

typedef struct {
  int count;
  int varindex;
} traverse_info;

static TriggerXGH *TriggerX_GetState(const cGH *GH) {
  return (TriggerXGH *)CCTK_GHExtension(GH, "TriggerX");
}

static void TriggerX_Traverse_Callback(int varindex, const char *optstring,
                                       void *arg) {
  traverse_info *info = (traverse_info *)arg;
  (void)optstring;
  info->count += 1;
  info->varindex = varindex;
}

static int TriggerX_ParseSingleVariable(const char *name, int *varindex) {
  traverse_info info;

  info.count = 0;
  info.varindex = -1;

  if (!name || CCTK_EQUALS(name, ""))
    return 0;

  if (!CCTK_TraverseString(name, TriggerX_Traverse_Callback, &info, CCTK_VAR))
    return 0;

  if (info.count != 1)
    return 0;

  *varindex = info.varindex;
  return 1;
}

static int TriggerX_GroupTypeFromVar(int varindex) {
  const int groupindex = CCTK_GroupIndexFromVarI(varindex);
  cGroup group;

  if (groupindex < 0)
    return -1;
  if (CCTK_GroupData(groupindex, &group))
    return -1;

  return group.grouptype;
}

static int TriggerX_CheckRelation(const char *relation, CCTK_REAL lhs,
                                  CCTK_REAL rhs) {
  if (CCTK_EQUALS(relation, ">"))
    return lhs > rhs;
  if (CCTK_EQUALS(relation, ">="))
    return lhs >= rhs;
  if (CCTK_EQUALS(relation, "<"))
    return lhs < rhs;
  if (CCTK_EQUALS(relation, "<="))
    return lhs <= rhs;
  if (CCTK_EQUALS(relation, "=="))
    return lhs == rhs;
  if (CCTK_EQUALS(relation, "!="))
    return lhs != rhs;

  CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
             "Unknown relation '%s'", relation);
  return 0;
}

static int TriggerX_UnpackReducedValue(int trigger, int vartype,
                                       const void *storage,
                                       CCTK_REAL *value) {
  switch (vartype) {
  case CCTK_VARIABLE_REAL:
    *value = *(const CCTK_REAL *)storage;
    return 1;
  case CCTK_VARIABLE_INT:
    *value = (CCTK_REAL)*(const CCTK_INT *)storage;
    return 1;
  default:
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Trigger %d uses reduction on unsupported variable type %d",
               trigger, vartype);
    return 0;
  }
}

static int TriggerX_CarpetXReductionSupported(const char *reduction) {
  return CCTK_EQUALS(reduction, "minimum") ||
         CCTK_EQUALS(reduction, "maximum") ||
         CCTK_EQUALS(reduction, "sum") ||
         CCTK_EQUALS(reduction, "sum_abs") ||
         CCTK_EQUALS(reduction, "sum_squared") ||
         CCTK_EQUALS(reduction, "average") ||
         CCTK_EQUALS(reduction, "standard_deviation") ||
         CCTK_EQUALS(reduction, "volume") ||
         CCTK_EQUALS(reduction, "norm1") ||
         CCTK_EQUALS(reduction, "norm2") ||
         CCTK_EQUALS(reduction, "norm_inf") ||
         CCTK_EQUALS(reduction, "maximum_abs");
}

static int TriggerX_CarpetXReductionValue(
    int trigger, const char *reduction, const CarpetX::reduction_CCTK_REAL &red,
    CCTK_REAL *value) {
  if (CCTK_EQUALS(reduction, "minimum")) {
    *value = red.min;
  } else if (CCTK_EQUALS(reduction, "maximum")) {
    *value = red.max;
  } else if (CCTK_EQUALS(reduction, "sum")) {
    *value = red.sum;
  } else if (CCTK_EQUALS(reduction, "sum_abs")) {
    *value = red.sumabs;
  } else if (CCTK_EQUALS(reduction, "sum_squared")) {
    *value = red.sum2;
  } else if (CCTK_EQUALS(reduction, "average")) {
    *value = red.avg();
  } else if (CCTK_EQUALS(reduction, "standard_deviation")) {
    *value = red.sdv();
  } else if (CCTK_EQUALS(reduction, "volume")) {
    *value = red.norm0();
  } else if (CCTK_EQUALS(reduction, "norm1")) {
    *value = red.norm1();
  } else if (CCTK_EQUALS(reduction, "norm2")) {
    *value = red.norm2();
  } else if (CCTK_EQUALS(reduction, "norm_inf") ||
             CCTK_EQUALS(reduction, "maximum_abs")) {
    *value = red.norm_inf();
  } else {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Trigger %d uses unsupported CarpetX reduction '%s'", trigger,
               reduction);
    return 0;
  }

  return 1;
}

static int TriggerX_ReduceVariable(const cGH *GH, int trigger, int varindex,
                                   const char *reduction, CCTK_REAL *value) {
  const int grouptype = TriggerX_GroupTypeFromVar(varindex);

  if (grouptype == CCTK_GF) {
    const int groupindex = CCTK_GroupIndexFromVarI(varindex);
    const int firstvarindex = CCTK_FirstVarIndexI(groupindex);
    const int vi = varindex - firstvarindex;

    if (groupindex < 0 || firstvarindex < 0 || vi < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Trigger %d could not inspect checked variable '%s'", trigger,
                 CCTK_FullVarName(varindex));
      return 0;
    }

    if (!TriggerX_CarpetXReductionSupported(reduction)) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Trigger %d uses unsupported reduction '%s' for grid function '%s'",
                 trigger, reduction, CCTK_FullVarName(varindex));
      return 0;
    }

    return TriggerX_CarpetXReductionValue(trigger, reduction,
                                          CarpetX::reduce(groupindex, vi, 0),
                                          value);
  }

  {
    const int vartype = CCTK_VarTypeI(varindex);
    const int reduction_handle = CCTK_ReductionHandle(reduction);

    if (reduction_handle < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Could not get reduction handle '%s' for trigger %d",
                 reduction, trigger);
      return 0;
    }

    if (vartype == CCTK_VARIABLE_REAL) {
      CCTK_REAL realval = 0.0;

      if (CCTK_Reduce(GH, -1, reduction_handle, 1, vartype, &realval, 1,
                      varindex)) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Reduction '%s' failed for trigger %d", reduction, trigger);
        return 0;
      }

      return TriggerX_UnpackReducedValue(trigger, vartype, &realval, value);
    }

    if (vartype == CCTK_VARIABLE_INT) {
      CCTK_INT intval = 0;

      if (CCTK_Reduce(GH, -1, reduction_handle, 1, vartype, &intval, 1,
                      varindex)) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Reduction '%s' failed for trigger %d", reduction, trigger);
        return 0;
      }

      return TriggerX_UnpackReducedValue(trigger, vartype, &intval, value);
    }

    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Reduction '%s' is unsupported for trigger %d variable type %d",
               reduction, trigger, vartype);
    return 0;
  }
}

static int TriggerX_ReadValue(const cGH *GH, TriggerXGH *my_GH, int trigger,
                              CCTK_REAL *value) {
  const int varindex = my_GH->checked_variable[trigger];
  int use_reduction = 0;

  if (!CCTK_EQUALS(my_GH->reduction[trigger], "")) {
    use_reduction = 1;
    if (varindex < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Trigger %d requests a reduction for a parameter; this is not "
                 "supported",
                 trigger);
      return 0;
    }
  }

  if (use_reduction) {
    return TriggerX_ReduceVariable(GH, trigger, varindex,
                                   my_GH->reduction[trigger], value);
  }

  if (varindex >= 0) {
    void *myVar = CCTK_VarDataPtrI(GH, 0, varindex);
    const int vartype = CCTK_VarTypeI(varindex);

    if (myVar == NULL) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Variable '%s' has no storage; skipping trigger %d",
                 CCTK_FullVarName(varindex), trigger);
      return 0;
    }

    if (vartype == CCTK_VARIABLE_REAL) {
      *value = *(const CCTK_REAL *)myVar;
      return 1;
    }
    if (vartype == CCTK_VARIABLE_INT) {
      *value = (CCTK_REAL)*(const CCTK_INT *)myVar;
      return 1;
    }

    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Variable '%s' is neither CCTK_REAL nor CCTK_INT",
               CCTK_FullVarName(varindex));
    return 0;
  }

  {
    int type;
    const void *tmp_value =
        CCTK_ParameterGet(my_GH->checked_parameter_name[trigger],
                          my_GH->checked_parameter_thorn[trigger], &type);

    if (!tmp_value) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Parameter '%s::%s' was not found",
                 my_GH->checked_parameter_thorn[trigger],
                 my_GH->checked_parameter_name[trigger]);
      return 0;
    }

    switch (type) {
    case PARAMETER_REAL:
      *value = *(const CCTK_REAL *)tmp_value;
      return 1;
    case PARAMETER_INT:
      *value = (CCTK_REAL)*(const CCTK_INT *)tmp_value;
      return 1;
    case PARAMETER_BOOLEAN:
      *value = (CCTK_REAL)*(const int *)tmp_value;
      return 1;
    default:
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Trigger %d uses unsupported parameter type %d", trigger,
                 type);
      return 0;
    }
  }
}

static int TriggerX_TriggerFulfilled(const cGH *GH, int trigger) {
  TriggerXGH *my_GH = TriggerX_GetState(GH);
  CCTK_REAL value = 0.0;
  const int have_value = TriggerX_ReadValue(GH, my_GH, trigger, &value);
  const int varindex = my_GH->checked_variable[trigger];
  const int ret =
      have_value &&
      TriggerX_CheckRelation(my_GH->relation[trigger], value,
                             my_GH->checked_value[trigger]);

  if (!my_GH->debug || !have_value)
    return ret;

  if (varindex >= 0) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Trigger %d %s for %s (%.17g %s %.17g)", trigger,
               ret ? "fulfilled" : "not fulfilled", CCTK_VarName(varindex),
               (double)value, my_GH->relation[trigger],
               (double)my_GH->checked_value[trigger]);
  } else {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Trigger %d %s for %s::%s (%.17g %s %.17g)", trigger,
               ret ? "fulfilled" : "not fulfilled",
               my_GH->checked_parameter_thorn[trigger],
               my_GH->checked_parameter_name[trigger], (double)value,
               my_GH->relation[trigger], (double)my_GH->checked_value[trigger]);
  }

  return ret;
}

static int TriggerX_SteerParameter(const char *parameter_name,
                                   const char *parameter_thorn,
                                   const char *parameter_value) {
  char *valstr =
      CCTK_ParameterValString(parameter_name, parameter_thorn);
  int ret;

  if (valstr && CCTK_EQUALS(valstr, parameter_value)) {
    free(valstr);
    return 0;
  }

  if (valstr)
    free(valstr);

  ret = CCTK_ParameterSet(parameter_name, parameter_thorn, parameter_value);

  if (ret == -1) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Parameter '%s::%s' is out of range",
               parameter_thorn, parameter_name);
  } else if (ret == -2) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Parameter '%s::%s' was not found",
               parameter_thorn, parameter_name);
  } else if (ret == -3) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Parameter '%s::%s' is not steerable",
               parameter_thorn, parameter_name);
  } else if (ret != 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not set parameter '%s::%s'",
               parameter_thorn, parameter_name);
  }

  return ret == 0 ? 1 : -1;
}

static int TriggerX_SteerScalar(const cGH *GH, TriggerXGH *my_GH, int trigger,
                                int scalar_index,
                                const char *scalar_value_string) {
  void *myVar;
  const int varindex = my_GH->steered_scalar[trigger];
  const int vartype = CCTK_VarTypeI(varindex);

  if (varindex < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Trigger %d has no scalar target", trigger);
    return -1;
  }

  myVar = CCTK_VarDataPtrI(GH, 0, varindex);
  if (myVar == NULL) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Variable '%s' has no storage; cannot steer trigger %d",
               CCTK_FullVarName(varindex), trigger);
    return -1;
  }

  if (vartype == CCTK_VARIABLE_REAL) {
    CCTK_REAL *real_var = (CCTK_REAL *)myVar;
    real_var[scalar_index] = atof(scalar_value_string);
    return 1;
  }

  if (vartype == CCTK_VARIABLE_INT) {
    CCTK_INT *int_var = (CCTK_INT *)myVar;
    int_var[scalar_index] = (CCTK_INT)atoi(scalar_value_string);
    return 1;
  }

  CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
             "Variable '%s' has unsupported type %d",
             CCTK_FullVarName(varindex), vartype);
  return -1;
}

extern "C" void TriggerX_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  TriggerXGH *my_GH = TriggerX_GetState(cctkGH);
  int i;

  trigger_cctk_iteration[0] = (CCTK_REAL)cctk_iteration;
  trigger_cctk_time[0] = (CCTK_REAL)cctk_time;

  if (!my_GH)
    CCTK_ERROR("Could not obtain TriggerX GH extension");

  for (i = 0; i < my_GH->number; ++i) {
    if (!my_GH->active[i])
      continue;

    my_GH->trigger_once[i] = Trigger_Once[i];

    if (my_GH->last_checked[i] >= cctk_iteration)
      continue;

    if (my_GH->trigger_count[i] > 0 && my_GH->trigger_once[i]) {
      my_GH->last_checked[i] = cctk_iteration;
      continue;
    }

    if (TriggerX_TriggerFulfilled(cctkGH, i)) {
      if (CCTK_EQUALS(my_GH->reaction[i], "steerparam")) {
        if (my_GH->debug)
          CCTK_VInfo(CCTK_THORNSTRING, "Trigger %d steering parameter", i);
        if (TriggerX_SteerParameter(Trigger_Steered_Parameter_Name[i],
                                    Trigger_Steered_Parameter_Thorn[i],
                                    Trigger_Steered_Parameter_Value[i]) > 0)
          my_GH->trigger_count[i] += 1;
      } else if (CCTK_EQUALS(my_GH->reaction[i], "steerscalar")) {
        if (my_GH->debug)
          CCTK_VInfo(CCTK_THORNSTRING, "Trigger %d steering scalar", i);
        if (TriggerX_SteerScalar(cctkGH, my_GH, i,
                                 Trigger_Steered_Scalar_Index[i],
                                 Trigger_Steered_Scalar_Value[i]) > 0)
          my_GH->trigger_count[i] += 1;
      }
    }

    my_GH->last_checked[i] = cctk_iteration;
  }
}

static void *TriggerX_SetupGH(tFleshConfig *config, int conv_level, cGH *GH) {
  DECLARE_CCTK_PARAMETERS;
  TriggerXGH *my_GH;
  int i;

  (void)config;
  (void)conv_level;

  my_GH = (TriggerXGH *)malloc(sizeof(TriggerXGH));

  my_GH->number = Trigger_Number;
  my_GH->debug = Trigger_Debug;
  my_GH->active = (int *)calloc(Trigger_Number, sizeof(int));
  my_GH->checked_variable = (int *)calloc(Trigger_Number, sizeof(int));
  my_GH->steered_scalar = (int *)calloc(Trigger_Number, sizeof(int));
  my_GH->last_checked = (int *)calloc(Trigger_Number, sizeof(int));
  my_GH->trigger_once = (int *)calloc(Trigger_Number, sizeof(int));
  my_GH->trigger_count = (int *)calloc(Trigger_Number, sizeof(int));
  my_GH->relation = (const char **)calloc(Trigger_Number, sizeof(const char *));
  my_GH->reduction = (const char **)calloc(Trigger_Number, sizeof(const char *));
  my_GH->checked_parameter_thorn =
      (const char **)calloc(Trigger_Number, sizeof(const char *));
  my_GH->checked_parameter_name =
      (const char **)calloc(Trigger_Number, sizeof(const char *));
  my_GH->reaction = (const char **)calloc(Trigger_Number, sizeof(const char *));
  my_GH->checked_value =
      (CCTK_REAL *)calloc(Trigger_Number, sizeof(CCTK_REAL));

  for (i = 0; i < Trigger_Number; ++i) {
    int cond_ok = 0;
    int target_ok = 0;

    my_GH->checked_variable[i] = -1;
    my_GH->steered_scalar[i] = -1;
    my_GH->last_checked[i] = -1;
    my_GH->trigger_once[i] = Trigger_Once[i];
    my_GH->trigger_count[i] = 0;
    my_GH->relation[i] = Trigger_Relation[i];
    my_GH->reduction[i] = Trigger_Reduction[i];
    my_GH->checked_value[i] = Trigger_Checked_Value[i];
    my_GH->reaction[i] = Trigger_Reaction[i];

    if (CCTK_EQUALS(Trigger_Checked_Variable[i], "param")) {
      if (CCTK_EQUALS(Trigger_Reduction[i], "")) {
        if (!CCTK_EQUALS(Trigger_Checked_Parameter_Name[i], "") &&
            !CCTK_EQUALS(Trigger_Checked_Parameter_Thorn[i], "") &&
            CCTK_ParameterGet(Trigger_Checked_Parameter_Name[i],
                              Trigger_Checked_Parameter_Thorn[i], NULL)) {
          my_GH->checked_parameter_name[i] = Trigger_Checked_Parameter_Name[i];
          my_GH->checked_parameter_thorn[i] = Trigger_Checked_Parameter_Thorn[i];
          cond_ok = 1;
        } else {
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Trigger %d refers to unknown parameter %s::%s", i,
                     Trigger_Checked_Parameter_Thorn[i],
                     Trigger_Checked_Parameter_Name[i]);
        }
      } else {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d cannot use a reduction on a parameter", i);
      }
    } else if (!CCTK_EQUALS(Trigger_Checked_Variable[i], "")) {
      if (TriggerX_ParseSingleVariable(Trigger_Checked_Variable[i],
                                       &my_GH->checked_variable[i])) {
        cond_ok = 1;
      } else {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d must reference exactly one checked variable: '%s'",
                   i, Trigger_Checked_Variable[i]);
      }
    }

    if (CCTK_EQUALS(Trigger_Reaction[i], "steerparam")) {
      if (CCTK_ParameterGet(Trigger_Steered_Parameter_Name[i],
                            Trigger_Steered_Parameter_Thorn[i], NULL)) {
        target_ok = 1;
      } else {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d refers to unknown steered parameter %s::%s", i,
                   Trigger_Steered_Parameter_Thorn[i],
                   Trigger_Steered_Parameter_Name[i]);
      }
    } else if (CCTK_EQUALS(Trigger_Reaction[i], "steerscalar")) {
      if (TriggerX_ParseSingleVariable(Trigger_Steered_Scalar[i],
                                       &my_GH->steered_scalar[i])) {
        const int grouptype = TriggerX_GroupTypeFromVar(my_GH->steered_scalar[i]);
        if (grouptype == CCTK_GF) {
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Trigger %d target '%s' is a grid function, not a scalar "
                     "or array",
                     i, Trigger_Steered_Scalar[i]);
        } else if (grouptype < 0) {
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Trigger %d could not inspect '%s'", i,
                     Trigger_Steered_Scalar[i]);
        } else {
          target_ok = 1;
        }
      } else {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d must reference exactly one steered scalar: '%s'",
                   i, Trigger_Steered_Scalar[i]);
      }
    }

    my_GH->active[i] = cond_ok && target_ok;
  }

  (void)GH;
  return my_GH;
}

extern "C" int TriggerX_Startup(void) {
  CCTK_RegisterGHExtensionSetupGH(CCTK_RegisterGHExtension("TriggerX"),
                                  TriggerX_SetupGH);
  return 0;
}

extern "C" void TriggerX_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  int i;

  for (i = 0; i < Trigger_Number; ++i) {
    if (CCTK_EQUALS(Trigger_Checked_Variable[i], "param")) {
      if (CCTK_EQUALS(Trigger_Checked_Parameter_Thorn[i], "") ||
          CCTK_EQUALS(Trigger_Checked_Parameter_Name[i], "")) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d uses 'param' but does not define the checked "
                   "parameter name and thorn",
                   i);
      }
      if (!CCTK_EQUALS(Trigger_Reduction[i], "")) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d cannot use a reduction on a parameter", i);
      }
    } else if (!CCTK_EQUALS(Trigger_Checked_Variable[i], "")) {
      int varindex = -1;

      if (!TriggerX_ParseSingleVariable(Trigger_Checked_Variable[i], &varindex)) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d must reference exactly one checked variable: '%s'",
                   i, Trigger_Checked_Variable[i]);
      } else if (CCTK_EQUALS(Trigger_Reduction[i], "") &&
                 TriggerX_GroupTypeFromVar(varindex) == CCTK_GF) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d checks grid function '%s' without a reduction; "
                   "this will use only the first local value",
                   i, Trigger_Checked_Variable[i]);
      } else if (!CCTK_EQUALS(Trigger_Reduction[i], "") &&
                 TriggerX_GroupTypeFromVar(varindex) == CCTK_GF &&
                 !TriggerX_CarpetXReductionSupported(Trigger_Reduction[i])) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d uses unsupported CarpetX reduction '%s' for '%s'",
                   i, Trigger_Reduction[i], Trigger_Checked_Variable[i]);
      }
    }

    if (CCTK_EQUALS(Trigger_Reaction[i], "steerparam")) {
      const cParamData *paramdata = CCTK_ParameterData(
          Trigger_Steered_Parameter_Name[i], Trigger_Steered_Parameter_Thorn[i]);
      if (!paramdata) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d steers unknown parameter %s::%s", i,
                   Trigger_Steered_Parameter_Thorn[i],
                   Trigger_Steered_Parameter_Name[i]);
      } else if (paramdata->steerable != CCTK_STEERABLE_ALWAYS) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d steers parameter %s::%s, but it is not "
                   "STEERABLE=ALWAYS",
                   i, Trigger_Steered_Parameter_Thorn[i],
                   Trigger_Steered_Parameter_Name[i]);
      }
    } else if (CCTK_EQUALS(Trigger_Reaction[i], "steerscalar")) {
      int varindex = -1;

      if (!TriggerX_ParseSingleVariable(Trigger_Steered_Scalar[i], &varindex)) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d must reference exactly one steered scalar: '%s'",
                   i, Trigger_Steered_Scalar[i]);
      } else if (TriggerX_GroupTypeFromVar(varindex) == CCTK_GF) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Trigger %d target '%s' is a grid function, not a scalar or "
                   "array",
                   i, Trigger_Steered_Scalar[i]);
      }
    }
  }
}
