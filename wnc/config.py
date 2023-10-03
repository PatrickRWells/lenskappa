from cosmap.config.models.sky import SkyCoord, Quantity
from cosmap.config.analysis import CosmapAnalysisParameters
from pydantic import Field, BaseModel

class LensParameters(BaseModel):
    lens_coordinate: SkyCoord
    source_redshift: float

class Filters(BaseModel):
    limiting_magnitude: float | list[float]
    filters: dict = {}
    class Config:
        arbitrary_types_allowed = True

class GeometryParameters(BaseModel):
    radius: Quantity
    inner_radius: Quantity


class Main(CosmapAnalysisParameters):
    n_samples: int = Field(1000, description="Number of samples to draw from the sky")
    lens_parameters: dict[str, LensParameters]
    geometry_parameters: GeometryParameters
    filters: Filters
