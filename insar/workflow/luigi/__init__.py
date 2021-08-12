import structlog
from insar.logs import TASK_LOGGER, STATUS_LOGGER, COMMON_PROCESSORS

# TBD: Should this just be part of insar.logs?
structlog.configure(processors=COMMON_PROCESSORS)
